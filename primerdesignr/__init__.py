"""
primerdesignr — Primer design and thermodynamic analysis

Backend: primer3-py (Tm, dimers) + seqfold (hairpin structure) + Mathews (hairpin second opinion)
"""

from .thermo import (
    analyze_primer,
    analyze_pair,
    cross_dimer_matrix,
    calc_tm,
    calc_hairpin,
    calc_homodimer,
    calc_heterodimer,
    PrimerReport,
    PairReport,
    TmResult,
    HairpinResult,
    DimerResult,
)

from .assembly import (
    design_gibson_overlap,
    check_golden_gate,
    check_overhang,
    recommend_assembly,
    GibsonOverlap,
    GoldenGateResult,
    OverhangCheck,
    AssemblyRecommendation,
)

__version__ = '0.3.0'
