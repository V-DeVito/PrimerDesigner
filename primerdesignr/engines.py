"""
Multi-engine orchestrator for primer thermodynamic analysis.

Engines:
    seqfold       — SantaLucia 2004 NN params (DNA), Zuker DP. Capped ~60 nt.
    mathews_hairpin — Custom Turner/Mathews 1999 RNA params (T→U mapping).
    primer3       — Primer3 hairpin Tm + dimer ΔG.
    viennarna     — OPTIONAL. ViennaRNA 2.7+. Handles long sequences efficiently
                    and provides MFE, suboptimal structures, and ensemble ΔG
                    using DNA Mathews 2004 or RNA Turner 2004 parameters.

ViennaRNA is optional; if it isn't installed the orchestrator returns
``vienna=None`` and callers must tolerate that key being absent.

Benchmark on 460 experimentally determined structures (Binet et al. 2023),
primer-length sequences (16-30 nt):
    seqfold (SantaLucia DNA)        — 36.2% exact match, F1 0.864
    ViennaRNA RNA (Turner 2004)     — 65.3% exact match, F1 0.939  ← best
    ViennaRNA DNA (Mathews 2004)    — 53.8% exact match, F1 0.891
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

# ── Optional ViennaRNA import ───────────────────────────────────────────
# ViennaRNA is GPL-licensed and may not be available in every install.
# Surface a clean boolean rather than letting the ImportError bubble up.
try:
    import RNA as _vienna  # noqa: N811 — ViennaRNA exposes itself as `RNA`
    HAS_VIENNARNA = True
except Exception:  # pragma: no cover — environment-dependent
    _vienna = None
    HAS_VIENNARNA = False


# ── Result types ────────────────────────────────────────────────────────

@dataclass
class ViennaResult:
    """ViennaRNA structure prediction for a single sequence."""
    dg: float                       # MFE ΔG (kcal/mol)
    dot_bracket: str                # MFE dot-bracket
    ensemble_dg: Optional[float] = None  # Partition-function ensemble ΔG (kcal/mol)
    suboptimal: List[Tuple[str, float]] = field(default_factory=list)
    param_set: str = 'dna_mathews_2004'  # 'dna_mathews_2004' | 'rna_turner_2004'


@dataclass
class ViennaCofoldResult:
    """ViennaRNA cofold (heterodimer) prediction."""
    dg: float                       # Cofold ΔG (kcal/mol)
    dot_bracket: str                # Combined struct, primers separated by '&'
    param_set: str = 'dna_mathews_2004'


# ── ViennaRNA wrappers ──────────────────────────────────────────────────

def _is_rna(seq: str) -> bool:
    """RNA if any U/u is present, else DNA."""
    return 'U' in seq.upper()


def _load_vienna_params(seq: str) -> str:
    """
    Load the appropriate parameter set for `seq`.

    NOTE: ViennaRNA's param loader is *global state* — every fold_compound
    created after this call uses the loaded params until something else
    overrides them. Always call this before constructing a fold_compound
    if you care about the parameter set.

    Returns the canonical name of the param set used.
    """
    if not HAS_VIENNARNA:
        raise RuntimeError("ViennaRNA is not installed")
    if _is_rna(seq):
        _vienna.params_load_RNA_Turner2004()
        return 'rna_turner_2004'
    _vienna.params_load_DNA_Mathews2004()
    return 'dna_mathews_2004'


def vienna_fold(
    seq: str,
    temp: float = 37.0,
    suboptimal_window_kcal: float = 0.0,
    compute_ensemble: bool = False,
) -> Optional[ViennaResult]:
    """
    Predict secondary structure using ViennaRNA.

    Args:
        seq: DNA or RNA sequence (T or U).
        temp: Temperature in °C.
        suboptimal_window_kcal: If >0, also compute suboptimal structures
            within this kcal/mol window of the MFE. Pass 5.0 for a typical
            "structures within 5 kcal/mol" view.
        compute_ensemble: If True, also compute the partition-function
            ensemble free energy. Slightly slower.

    Returns ``None`` if ViennaRNA is not available.
    """
    if not HAS_VIENNARNA:
        return None

    seq = seq.upper()
    param_set = _load_vienna_params(seq)

    md = _vienna.md()
    md.temperature = temp
    fc = _vienna.fold_compound(seq, md)

    # MFE
    mfe_struct, mfe_energy = fc.mfe()

    ensemble_dg: Optional[float] = None
    if compute_ensemble:
        # pf() returns (centroid_struct, ensemble_free_energy_kcal/mol)
        try:
            _, ensemble_dg = fc.pf()
        except Exception:  # pragma: no cover — defensive on edge cases
            ensemble_dg = None

    suboptimal: List[Tuple[str, float]] = []
    if suboptimal_window_kcal and suboptimal_window_kcal > 0:
        # ViennaRNA `subopt(delta)` takes delta in 0.01 kcal/mol increments,
        # but the returned `s.energy` is already kcal/mol. Do NOT divide by 100.
        delta = int(round(suboptimal_window_kcal * 100))
        try:
            results = fc.subopt(delta)
            for s in results:
                if s.structure is None:
                    continue
                suboptimal.append((s.structure, float(s.energy)))
        except Exception:  # pragma: no cover
            suboptimal = []

    return ViennaResult(
        dg=round(float(mfe_energy), 2),
        dot_bracket=mfe_struct,
        ensemble_dg=round(ensemble_dg, 2) if ensemble_dg is not None else None,
        suboptimal=suboptimal,
        param_set=param_set,
    )


def vienna_cofold(
    seq1: str,
    seq2: str,
    temp: float = 37.0,
) -> Optional[ViennaCofoldResult]:
    """
    Predict heterodimer (two-sequence cofold) using ViennaRNA.

    Returns ``None`` if ViennaRNA is not available.
    """
    if not HAS_VIENNARNA:
        return None

    s1, s2 = seq1.upper(), seq2.upper()
    param_set = _load_vienna_params(s1 + s2)

    md = _vienna.md()
    md.temperature = temp
    combined = f"{s1}&{s2}"
    fc = _vienna.fold_compound(combined, md)
    structure, energy = fc.mfe_dimer()
    return ViennaCofoldResult(
        dg=round(float(energy), 2),
        dot_bracket=structure,
        param_set=param_set,
    )


__all__ = [
    'HAS_VIENNARNA',
    'ViennaResult',
    'ViennaCofoldResult',
    'vienna_fold',
    'vienna_cofold',
]
