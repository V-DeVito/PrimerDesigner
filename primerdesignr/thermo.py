"""
primerdesignr.thermo — Thermodynamic analysis engine

Backend architecture:
    primer3-py  →  Tm, homodimer ΔG, heterodimer ΔG  (SantaLucia 2004)
    seqfold     →  Hairpin MFE + dot-bracket structure (Zuker DP, SantaLucia 2004)
    mathews     →  Hairpin second opinion              (Turner/Mathews RNA params)

The Mathews hairpin check catches AT-closing stems that SantaLucia
over-penalizes — see Binet et al. BMC Bioinformatics (2023) 24:422.

primer3-py is GPL v2, contained server-side. seqfold is MIT.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple
import primer3
from seqfold import dg as seqfold_dg, fold as seqfold_fold, dot_bracket

from . import mathews_hairpin


# =============================================================================
# RESULT TYPES
# =============================================================================

@dataclass
class TmResult:
    """Melting temperature calculation result."""
    tm: float           # °C
    method: str = 'santalucia_nn'


@dataclass
class DimerResult:
    """Homodimer or heterodimer result."""
    dg: float           # ΔG in kcal/mol
    dg_threshold: float = -9.0  # Industry standard cutoff
    ascii_structure: str = ''

    @property
    def is_problematic(self) -> bool:
        return self.dg < self.dg_threshold


@dataclass
class HairpinResult:
    """Hairpin analysis with dual-engine comparison."""
    dg_santalucia: float        # From seqfold (SantaLucia 2004)
    dg_mathews: float           # From Mathews/Turner RNA params
    dot_bracket: str            # Structure notation from seqfold
    dg_threshold: float = -3.0  # Hairpin concern threshold

    @property
    def engines_disagree(self) -> bool:
        """True when SantaLucia says 'no hairpin' but Mathews says 'hairpin'."""
        return self.dg_santalucia > self.dg_threshold and self.dg_mathews < self.dg_threshold

    @property
    def is_problematic(self) -> bool:
        """Conservative: flag if either engine says problematic."""
        return self.dg_santalucia < self.dg_threshold or self.dg_mathews < self.dg_threshold

    @property
    def worst_dg(self) -> float:
        return min(self.dg_santalucia, self.dg_mathews)


@dataclass
class PrimerReport:
    """Complete analysis of a single primer."""
    sequence: str
    length: int
    gc_content: float
    tm: TmResult
    hairpin: HairpinResult
    homodimer: DimerResult
    warnings: List[str] = field(default_factory=list)


@dataclass
class PairReport:
    """Analysis of a forward/reverse primer pair."""
    forward: PrimerReport
    reverse: PrimerReport
    heterodimer: DimerResult
    tm_difference: float
    warnings: List[str] = field(default_factory=list)


# =============================================================================
# CORE CALCULATIONS
# =============================================================================

def calc_tm(
    seq: str,
    mv_conc: float = 50.0,     # Monovalent cation mM (Na+/K+)
    dv_conc: float = 0.0,      # Divalent cation mM (Mg2+) — IDT default is 0
    dntp_conc: float = 0.0,    # dNTP mM — IDT default is 0
    dna_conc: float = 250.0,   # Primer nM
) -> TmResult:
    """
    Calculate Tm using SantaLucia nearest-neighbor method.

    Default conditions match IDT OligoAnalyzer: 50mM Na+, no Mg, no dNTPs.
    For realistic PCR conditions, pass dv_conc=1.5, dntp_conc=0.2.

    Wraps primer3.calc_tm with SantaLucia salt correction.
    """
    tm = primer3.calc_tm(
        seq,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        tm_method='santalucia',
        salt_corrections_method='santalucia',
    )
    return TmResult(tm=round(tm, 1))


def calc_hairpin(seq: str, temp: float = 37.0) -> HairpinResult:
    """
    Dual-engine hairpin prediction.

    Returns both SantaLucia (seqfold) and Mathews (custom) ΔG values.
    When they disagree on AT-closing stems, the UI flags the discrepancy.
    """
    seq = seq.upper()

    # Engine 1: seqfold (SantaLucia 2004, Zuker DP)
    dg_sl = seqfold_dg(seq, temp=temp)
    structs = seqfold_fold(seq)
    db = dot_bracket(seq, structs)

    # Engine 2: Mathews/Turner RNA params (catches AT-closing stems)
    dg_mw = mathews_hairpin.calc_hairpin_dg(seq, temp=temp)

    return HairpinResult(
        dg_santalucia=round(dg_sl, 2),
        dg_mathews=round(dg_mw, 2),
        dot_bracket=db,
    )


def calc_homodimer(seq: str, temp: float = 37.0) -> DimerResult:
    """
    Homodimer (self-dimer) ΔG via primer3-py.

    Two copies of the same primer hybridizing to each other.
    """
    result = primer3.calc_homodimer(
        seq,
        mv_conc=50.0,
        dv_conc=1.5,
        temp_c=temp,
    )
    return DimerResult(
        dg=round(result.dg / 1000, 2),  # primer3 returns cal/mol
        ascii_structure=result.ascii_structure if hasattr(result, 'ascii_structure') else '',
    )


def calc_heterodimer(seq1: str, seq2: str, temp: float = 37.0) -> DimerResult:
    """
    Heterodimer (cross-dimer) ΔG via primer3-py.

    Forward and reverse primers hybridizing to each other.
    """
    result = primer3.calc_heterodimer(
        seq1, seq2,
        mv_conc=50.0,
        dv_conc=1.5,
        temp_c=temp,
    )
    return DimerResult(
        dg=round(result.dg / 1000, 2),
        ascii_structure=result.ascii_structure if hasattr(result, 'ascii_structure') else '',
    )


# =============================================================================
# HIGH-LEVEL ANALYSIS
# =============================================================================

def _gc_content(seq: str) -> float:
    seq = seq.upper()
    return sum(1 for b in seq if b in 'GC') / len(seq) if seq else 0.0


def _generate_warnings(seq: str, tm: TmResult, hairpin: HairpinResult,
                       homodimer: DimerResult) -> List[str]:
    """Generate actionable warnings for a primer."""
    warnings = []
    seq = seq.upper()
    n = len(seq)
    gc = _gc_content(seq)

    # Length
    if n < 18:
        warnings.append(f"Short ({n}nt) — recommend 18–25")
    elif n > 30:
        warnings.append(f"Long ({n}nt) — may reduce specificity")

    # GC content
    if gc < 0.40:
        warnings.append(f"Low GC ({gc:.0%}) — recommend 40–60%")
    elif gc > 0.60:
        warnings.append(f"High GC ({gc:.0%}) — recommend 40–60%")

    # Tm range
    if tm.tm < 52:
        warnings.append(f"Low Tm ({tm.tm}°C) — may not anneal at standard temps")
    elif tm.tm > 68:
        warnings.append(f"High Tm ({tm.tm}°C) — consider shortening primer")

    # Hairpin
    if hairpin.is_problematic:
        warnings.append(
            f"Hairpin concern (SL: {hairpin.dg_santalucia}, MW: {hairpin.dg_mathews} kcal/mol)"
        )
    if hairpin.engines_disagree:
        warnings.append(
            "⚠ SantaLucia/Mathews disagree on hairpin — likely AT-closing stem"
        )

    # Self-dimer
    if homodimer.is_problematic:
        warnings.append(f"Self-dimer ΔG = {homodimer.dg} kcal/mol (threshold: -9)")

    # 3' end: GC clamp check (want 1-2 G/C in last 5 bases, not 0, not 5)
    three_prime = seq[-5:]
    gc_3p = sum(1 for b in three_prime if b in 'GC')
    if gc_3p == 0:
        warnings.append("No GC in last 5 bases — weak 3' anchoring")
    elif gc_3p >= 4:
        warnings.append("Heavy GC at 3' end — risk of mispriming")

    # Homopolymer runs
    for base in 'ATGC':
        if base * 4 in seq:
            warnings.append(f"Run of 4+ {base}s")

    return warnings


def analyze_primer(
    seq: str,
    mv_conc: float = 50.0,
    dv_conc: float = 0.0,
    dntp_conc: float = 0.0,
    dna_conc: float = 250.0,
) -> PrimerReport:
    """Complete thermodynamic analysis of a single primer."""
    seq = seq.upper().strip()

    tm = calc_tm(seq, mv_conc, dv_conc, dntp_conc, dna_conc)
    hairpin = calc_hairpin(seq)
    homodimer = calc_homodimer(seq)
    warnings = _generate_warnings(seq, tm, hairpin, homodimer)

    return PrimerReport(
        sequence=seq,
        length=len(seq),
        gc_content=round(_gc_content(seq), 3),
        tm=tm,
        hairpin=hairpin,
        homodimer=homodimer,
        warnings=warnings,
    )


def analyze_pair(fwd: str, rev: str, **kwargs) -> PairReport:
    """Analyze a primer pair including cross-dimer check."""
    fwd_report = analyze_primer(fwd, **kwargs)
    rev_report = analyze_primer(rev, **kwargs)
    heterodimer = calc_heterodimer(fwd, rev)

    tm_diff = round(abs(fwd_report.tm.tm - rev_report.tm.tm), 1)

    warnings = []
    if tm_diff > 5:
        warnings.append(f"Tm difference {tm_diff}°C — recommend <5°C")
    if heterodimer.is_problematic:
        warnings.append(f"Cross-dimer ΔG = {heterodimer.dg} kcal/mol")

    return PairReport(
        forward=fwd_report,
        reverse=rev_report,
        heterodimer=heterodimer,
        tm_difference=tm_diff,
        warnings=warnings,
    )


def cross_dimer_matrix(primers: dict[str, str]) -> dict[tuple[str, str], DimerResult]:
    """
    Compute all pairwise heterodimer ΔG values.

    Args:
        primers: dict of {name: sequence}

    Returns:
        dict of {(name1, name2): DimerResult} for all unique pairs
    """
    names = list(primers.keys())
    results = {}

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            result = calc_heterodimer(primers[n1], primers[n2])
            results[(n1, n2)] = result

    return results
