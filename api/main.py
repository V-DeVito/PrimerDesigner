"""
primerdesignr API — FastAPI backend

Endpoints:
    POST /analyze          — analyze one or more primers
    POST /pair             — analyze a primer pair
    POST /cross-dimers     — full pairwise cross-dimer matrix
    POST /golden-gate      — validate Golden Gate overhangs
    GET  /health           — healthcheck
"""

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field, field_validator
from typing import Optional
import re
import sys
import os

# Add the parent directory so we can import primerdesignr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from primerdesignr.thermo import (
    analyze_primer, analyze_pair, cross_dimer_matrix,
)
from primerdesignr.assembly import check_golden_gate

app = FastAPI(
    title="primerdesignr",
    version="0.1.0",
    description="Primer design and thermodynamic analysis for synthetic biology",
)

# CORS — allow the SvelteKit frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",       # SvelteKit dev
        "http://localhost:4173",       # SvelteKit preview
        "https://primerdesignr.com",   # production
        "https://*.vercel.app",        # Vercel preview deploys
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Validation ──────────────────────────────────────────────

DNA_RE = re.compile(r'^[ATGCatgc]+$')
MAX_SEQ_LEN = 200
MAX_PRIMERS = 24  # 24 primers = 276 pairwise comparisons, still fast


def validate_dna(seq: str, field_name: str = "sequence") -> str:
    seq = seq.strip().upper()
    if not seq:
        raise HTTPException(400, f"{field_name} is empty")
    if not DNA_RE.match(seq):
        raise HTTPException(400, f"{field_name} contains non-ATGC characters")
    if len(seq) > MAX_SEQ_LEN:
        raise HTTPException(400, f"{field_name} exceeds {MAX_SEQ_LEN}nt limit")
    if len(seq) < 10:
        raise HTTPException(400, f"{field_name} is too short (<10nt)")
    return seq


# ── Request / Response Models ───────────────────────────────

class Conditions(BaseModel):
    na_mm: float = Field(50.0, ge=0, le=1000, description="Na+ concentration (mM)")
    mg_mm: float = Field(0.0, ge=0, le=100, description="Mg2+ concentration (mM)")
    dna_nm: float = Field(250.0, ge=1, le=10000, description="Primer concentration (nM)")


class AnalyzeRequest(BaseModel):
    primers: dict[str, str] = Field(
        ...,
        description="Map of primer name → sequence",
        min_length=1,
    )
    conditions: Conditions = Conditions()
    include_cross_dimers: bool = Field(True, description="Compute pairwise cross-dimers")


class PairRequest(BaseModel):
    forward: str
    reverse: str
    forward_name: str = "Forward"
    reverse_name: str = "Reverse"
    conditions: Conditions = Conditions()


class GoldenGateRequest(BaseModel):
    overhangs: list[str] = Field(..., min_length=2, max_length=20)
    enzyme: str = "BsaI"
    sequences: Optional[dict[str, str]] = None  # For internal site scanning


class HairpinResponse(BaseModel):
    dg_santalucia: float
    dg_mathews: float
    dot_bracket: str
    engines_disagree: bool
    is_problematic: bool


class DimerResponse(BaseModel):
    dg: float
    is_problematic: bool


class PrimerResponse(BaseModel):
    sequence: str
    length: int
    gc_content: float
    tm: float
    hairpin: HairpinResponse
    homodimer: DimerResponse
    warnings: list[str]


class CrossDimerEntry(BaseModel):
    primer_a: str
    primer_b: str
    dg: float
    is_problematic: bool


class AnalyzeResponse(BaseModel):
    primers: dict[str, PrimerResponse]
    cross_dimers: list[CrossDimerEntry]
    tm_spread: float
    summary: str


# ── Helpers ─────────────────────────────────────────────────

def _primer_to_response(report) -> PrimerResponse:
    return PrimerResponse(
        sequence=report.sequence,
        length=report.length,
        gc_content=report.gc_content,
        tm=report.tm.tm,
        hairpin=HairpinResponse(
            dg_santalucia=report.hairpin.dg_santalucia,
            dg_mathews=report.hairpin.dg_mathews,
            dot_bracket=report.hairpin.dot_bracket,
            engines_disagree=report.hairpin.engines_disagree,
            is_problematic=report.hairpin.is_problematic,
        ),
        homodimer=DimerResponse(
            dg=report.homodimer.dg,
            is_problematic=report.homodimer.is_problematic,
        ),
        warnings=report.warnings,
    )


# ── Endpoints ───────────────────────────────────────────────

@app.get("/health")
def health():
    return {"status": "ok", "version": "0.1.0"}


@app.post("/analyze", response_model=AnalyzeResponse)
def analyze(req: AnalyzeRequest):
    if len(req.primers) > MAX_PRIMERS:
        raise HTTPException(400, f"Maximum {MAX_PRIMERS} primers at once")

    # Validate sequences
    clean = {}
    for name, seq in req.primers.items():
        clean[name] = validate_dna(seq, name)

    # Analyze each primer
    reports = {}
    for name, seq in clean.items():
        reports[name] = analyze_primer(
            seq,
            mv_conc=req.conditions.na_mm,
            dv_conc=req.conditions.mg_mm,
            dna_conc=req.conditions.dna_nm,
        )

    # Cross-dimer matrix
    cross_dimers = []
    if req.include_cross_dimers and len(clean) > 1:
        matrix = cross_dimer_matrix(clean)
        for (n1, n2), result in matrix.items():
            cross_dimers.append(CrossDimerEntry(
                primer_a=n1,
                primer_b=n2,
                dg=result.dg,
                is_problematic=result.is_problematic,
            ))

    # Tm spread
    tms = [r.tm.tm for r in reports.values()]
    tm_spread = round(max(tms) - min(tms), 1) if len(tms) > 1 else 0.0

    # Summary
    problems = []
    problematic_dimers = [cd for cd in cross_dimers if cd.is_problematic]
    disagreements = [n for n, r in reports.items() if r.hairpin.engines_disagree]

    if problematic_dimers:
        problems.append(f"{len(problematic_dimers)} problematic cross-dimer(s)")
    if disagreements:
        problems.append(f"{len(disagreements)} hairpin disagreement(s)")
    if tm_spread > 5:
        problems.append(f"Tm spread {tm_spread}°C exceeds 5°C")

    summary = "; ".join(problems) if problems else "All clear"

    return AnalyzeResponse(
        primers={name: _primer_to_response(r) for name, r in reports.items()},
        cross_dimers=cross_dimers,
        tm_spread=tm_spread,
        summary=summary,
    )


@app.post("/pair")
def pair(req: PairRequest):
    fwd = validate_dna(req.forward, "forward")
    rev = validate_dna(req.reverse, "reverse")

    result = analyze_pair(
        fwd, rev,
        mv_conc=req.conditions.na_mm,
        dv_conc=req.conditions.mg_mm,
        dna_conc=req.conditions.dna_nm,
    )

    return {
        "forward": _primer_to_response(result.forward),
        "reverse": _primer_to_response(result.reverse),
        "heterodimer": DimerResponse(
            dg=result.heterodimer.dg,
            is_problematic=result.heterodimer.is_problematic,
        ),
        "tm_difference": result.tm_difference,
        "warnings": result.warnings,
    }


@app.post("/golden-gate")
def golden_gate(req: GoldenGateRequest):
    # Validate overhangs
    for oh in req.overhangs:
        oh = oh.strip().upper()
        if len(oh) != 4 or not DNA_RE.match(oh):
            raise HTTPException(400, f"Invalid overhang: {oh} (must be 4bp DNA)")

    result = check_golden_gate(
        req.overhangs,
        sequences=req.sequences,
        enzyme=req.enzyme,
    )

    return {
        "overhangs": [
            {
                "overhang": oh.overhang,
                "reverse_complement": oh.reverse_complement,
                "gc_content": oh.gc_content,
                "is_palindromic": oh.is_palindromic,
                "is_in_fidelity_set": oh.is_in_fidelity_set,
                "warnings": oh.warnings,
            }
            for oh in result.overhangs
        ],
        "all_unique": result.all_unique,
        "enzyme": result.enzyme,
        "internal_sites": [
            {"sequence": name, "position": pos}
            for name, pos in result.internal_sites
        ],
        "warnings": result.warnings,
    }
