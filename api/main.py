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
from typing import Literal, Optional
import re
import sys
import os

# Add the parent directory so we can import primerdesignr
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from primerdesignr.thermo import (
    analyze_primer, analyze_pair, cross_dimer_matrix,
)
from primerdesignr.assembly import check_golden_gate
from primerdesignr.design import design_pcr_primers

app = FastAPI(
    title="primerdesignr",
    version="0.1.0",
    description="Primer design and thermodynamic analysis for synthetic biology",
)

def _cors_origins() -> list[str]:
    raw = os.getenv("CORS_ORIGINS", "")
    if raw:
        return [origin.strip() for origin in raw.split(",") if origin.strip()]

    return [
        "http://localhost:5173",
        "http://127.0.0.1:5173",
        "http://localhost:4173",
        "http://127.0.0.1:4173",
    ]


# CORS is intentionally an allowlist. Set CORS_ORIGINS in production.
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins(),
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Validation ──────────────────────────────────────────────

DNA_RE = re.compile(r'^[ATGCatgc]+$')
MAX_SEQ_LEN = 200
MAX_TEMPLATE_LEN = 20000
MAX_PRIMERS = 24  # 24 primers = 276 pairwise comparisons, still fast


def clean_dna_text(seq: str) -> str:
    """Accept plain DNA or FASTA-like pasted text and strip headers/whitespace."""
    lines = []
    for line in seq.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(">"):
            continue
        lines.append(stripped)
    return re.sub(r"\s+", "", "".join(lines)).upper()


def validate_dna(
    seq: str,
    field_name: str = "sequence",
    min_len: int = 10,
    max_len: int = MAX_SEQ_LEN,
) -> str:
    seq = clean_dna_text(seq)
    if not seq:
        raise HTTPException(400, f"{field_name} is empty")
    if not DNA_RE.match(seq):
        raise HTTPException(400, f"{field_name} contains non-ATGC characters")
    if len(seq) > max_len:
        raise HTTPException(400, f"{field_name} exceeds {max_len}nt limit")
    if len(seq) < min_len:
        raise HTTPException(400, f"{field_name} is too short (<{min_len}nt)")
    return seq


# ── Request / Response Models ───────────────────────────────

class Conditions(BaseModel):
    na_mm: float = Field(50.0, ge=0, le=1000, description="Na+ concentration (mM)")
    mg_mm: float = Field(0.0, ge=0, le=100, description="Mg2+ concentration (mM)")
    dntp_mm: float = Field(0.0, ge=0, le=20, description="dNTP concentration (mM)")
    dna_nm: float = Field(250.0, ge=1, le=10000, description="Primer concentration (nM)")


class AnalyzeRequest(BaseModel):
    primers: dict[str, str] = Field(
        ...,
        description="Map of primer name → sequence",
        min_length=1,
    )
    conditions: Conditions = Field(default_factory=Conditions)
    include_cross_dimers: bool = Field(True, description="Compute pairwise cross-dimers")

    @field_validator("primers")
    @classmethod
    def validate_primer_names(cls, value: dict[str, str]) -> dict[str, str]:
        cleaned = {}
        for name, seq in value.items():
            clean_name = str(name).strip()
            if not clean_name:
                raise ValueError("Primer names cannot be blank")
            if len(clean_name) > 80:
                raise ValueError(f"Primer name is too long: {clean_name[:40]}")
            cleaned[clean_name] = seq
        return cleaned


class PairRequest(BaseModel):
    forward: str
    reverse: str
    forward_name: str = "Forward"
    reverse_name: str = "Reverse"
    conditions: Conditions = Field(default_factory=Conditions)


class DesignRequest(BaseModel):
    template: str = Field(..., min_length=20, description="Template DNA or FASTA text")
    design_mode: Literal["exact", "targeted", "amplicon"] = Field("exact", description="Exact bounds, required-region, or exploratory internal design")
    product_min: int = Field(120, ge=40, le=5000)
    product_max: int = Field(500, ge=40, le=5000)
    primer_count: int = Field(5, ge=1, le=20)
    target_start: Optional[int] = Field(None, ge=0)
    target_length: Optional[int] = Field(None, ge=1, le=5000)
    conditions: Conditions = Field(default_factory=Conditions)

    @field_validator("product_max")
    @classmethod
    def validate_product_max(cls, value: int, info):
        product_min = info.data.get("product_min")
        if product_min is not None and value < product_min:
            raise ValueError("product_max must be greater than or equal to product_min")
        return value


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


class PrimerCoordinatesResponse(BaseModel):
    start: int
    end: int
    length: int
    strand: str


class DesignCandidateResponse(BaseModel):
    rank: int
    forward: PrimerResponse
    reverse: PrimerResponse
    forward_coords: PrimerCoordinatesResponse
    reverse_coords: PrimerCoordinatesResponse
    product_size: int
    heterodimer: DimerResponse
    tm_difference: float
    primer3_pair_penalty: float
    primer3_pair_compl_any: Optional[float]
    primer3_pair_compl_end: Optional[float]
    explanations: list[str]
    warnings: list[str]


class DesignResponse(BaseModel):
    template_length: int
    product_size_range: tuple[int, int]
    target_region: Optional[tuple[int, int]]
    candidates: list[DesignCandidateResponse]
    primer3_explain: dict[str, str]
    warnings: list[str]


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


def _coords_to_response(coords) -> PrimerCoordinatesResponse:
    return PrimerCoordinatesResponse(
        start=coords.start,
        end=coords.end,
        length=coords.length,
        strand=coords.strand,
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
            dntp_conc=req.conditions.dntp_mm,
            dna_conc=req.conditions.dna_nm,
        )

    # Cross-dimer matrix
    cross_dimers = []
    if req.include_cross_dimers and len(clean) > 1:
        matrix = cross_dimer_matrix(
            clean,
            mv_conc=req.conditions.na_mm,
            dv_conc=req.conditions.mg_mm,
            dntp_conc=req.conditions.dntp_mm,
            dna_conc=req.conditions.dna_nm,
        )
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
        dntp_conc=req.conditions.dntp_mm,
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


@app.post("/design", response_model=DesignResponse)
def design(req: DesignRequest):
    template = validate_dna(
        req.template,
        "template",
        min_len=max(40, req.product_min) if req.design_mode in {"targeted", "amplicon"} else 40,
        max_len=MAX_TEMPLATE_LEN,
    )
    if req.design_mode in {"targeted", "amplicon"} and req.product_max > len(template):
        raise HTTPException(400, "product_max cannot exceed template length")
    if req.design_mode == "targeted" and (req.target_start is None or req.target_length is None):
        raise HTTPException(400, "target_start and target_length are required for required-region design")
    if req.target_start is None and req.target_length is not None:
        raise HTTPException(400, "target_start is required when target_length is set")
    if req.target_start is not None and req.target_length is None:
        raise HTTPException(400, "target_length is required when target_start is set")
    if req.target_start is not None and req.target_length is not None:
        if req.target_start + req.target_length > len(template):
            raise HTTPException(400, "target region extends beyond template length")

    result = design_pcr_primers(
        template,
        product_min=req.product_min,
        product_max=req.product_max,
        primer_count=req.primer_count,
        target_start=req.target_start,
        target_length=req.target_length,
        mv_conc=req.conditions.na_mm,
        dv_conc=req.conditions.mg_mm,
        dntp_conc=req.conditions.dntp_mm,
        dna_conc=req.conditions.dna_nm,
        design_mode=req.design_mode,
    )

    return DesignResponse(
        template_length=result.template_length,
        product_size_range=result.product_size_range,
        target_region=result.target_region,
        candidates=[
            DesignCandidateResponse(
                rank=c.rank,
                forward=_primer_to_response(c.pair.forward),
                reverse=_primer_to_response(c.pair.reverse),
                forward_coords=_coords_to_response(c.forward_coords),
                reverse_coords=_coords_to_response(c.reverse_coords),
                product_size=c.product_size,
                heterodimer=DimerResponse(
                    dg=c.heterodimer.dg,
                    is_problematic=c.heterodimer.is_problematic,
                ),
                tm_difference=c.pair.tm_difference,
                primer3_pair_penalty=c.primer3_pair_penalty,
                primer3_pair_compl_any=c.primer3_pair_compl_any,
                primer3_pair_compl_end=c.primer3_pair_compl_end,
                explanations=c.explanations,
                warnings=c.warnings,
            )
            for c in result.candidates
        ],
        primer3_explain=result.primer3_explain,
        warnings=result.warnings,
    )


@app.post("/golden-gate")
def golden_gate(req: GoldenGateRequest):
    # Validate overhangs
    for oh in req.overhangs:
        oh = oh.strip().upper()
        if len(oh) != 4 or not DNA_RE.match(oh):
            raise HTTPException(400, f"Invalid overhang: {oh} (must be 4bp DNA)")

    result = check_golden_gate(
        req.overhangs,
        sequences={
            name: validate_dna(seq, name, min_len=1, max_len=MAX_TEMPLATE_LEN)
            for name, seq in req.sequences.items()
        } if req.sequences else None,
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
