"""
primerdesignr.assembly — Assembly-specific primer design logic

Gibson Assembly:
    - Overlap regions need Tm 48–65°C (NEB recommends ~50°C for 20–40bp overlaps)
    - Binding region Tm should be 55–72°C
    - Overlap length: 15–40bp, recommend 20–30bp

Golden Gate:
    - 4bp overhangs must be unique and non-palindromic
    - Checks for ligation fidelity using NEB's fidelity data
    - BsaI (GGTCTC) / BsmBI (CGTCTC) recognition site screening

References:
    - Gibson DG et al. (2009) Nature Methods 6:343-345
    - Potapov V et al. (2018) ACS Synth Biol 7:2665-2674 (Golden Gate fidelity)
    - NEB Golden Gate Assembly Tool guidelines
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
import primer3

from .thermo import calc_tm, TmResult


# =============================================================================
# GIBSON ASSEMBLY
# =============================================================================

@dataclass
class GibsonOverlap:
    """Analysis of a Gibson assembly overlap region."""
    overlap_seq: str
    overlap_tm: TmResult
    binding_seq: str
    binding_tm: TmResult
    full_primer: str
    overlap_length: int
    warnings: List[str] = field(default_factory=list)


def design_gibson_overlap(
    upstream_3prime: str,
    downstream_5prime: str,
    binding_region: str,
    overlap_target_tm: float = 50.0,
    binding_target_tm: float = 60.0,
    min_overlap: int = 15,
    max_overlap: int = 40,
    mv_conc: float = 50.0,
    dv_conc: float = 0.0,
    dna_conc: float = 250.0,
) -> GibsonOverlap:
    """
    Design a Gibson assembly primer with optimized overlap Tm.

    The primer has two parts:
        5'-[OVERLAP (anneals to adjacent fragment)]--[BINDING (anneals to template)]--3'

    The overlap Tm is calculated for the overlap region annealing to the
    adjacent fragment during the isothermal assembly step (50°C).
    The binding Tm is for the initial PCR amplification step.

    Args:
        upstream_3prime: 3' end of the upstream fragment (overlap source)
        downstream_5prime: 5' end of the downstream fragment (overlap source)
        binding_region: Template region the primer binds to for PCR
        overlap_target_tm: Target Tm for overlap (default 50°C for Gibson)
        binding_target_tm: Target Tm for binding region
        min_overlap: Minimum overlap length
        max_overlap: Maximum overlap length
    """
    warnings = []

    # Build overlap from the junction of two fragments
    # For a forward primer at a junction: overlap = rc(upstream_3prime)
    # This is context-dependent; here we just take the provided overlap sequence
    overlap_source = upstream_3prime[-max_overlap:] + downstream_5prime[:max_overlap]

    # Find optimal overlap length by scanning for target Tm
    best_overlap = None
    best_tm_diff = float('inf')

    for length in range(min_overlap, min(max_overlap + 1, len(upstream_3prime) + 1)):
        candidate = upstream_3prime[-length:]
        tm = calc_tm(candidate, mv_conc=mv_conc, dv_conc=dv_conc, dna_conc=dna_conc)
        diff = abs(tm.tm - overlap_target_tm)

        if diff < best_tm_diff:
            best_tm_diff = diff
            best_overlap = candidate
            best_overlap_tm = tm

    if best_overlap is None:
        best_overlap = upstream_3prime[-min_overlap:]
        best_overlap_tm = calc_tm(best_overlap, mv_conc=mv_conc, dv_conc=dv_conc, dna_conc=dna_conc)

    # Calculate binding region Tm
    binding_tm = calc_tm(binding_region, mv_conc=mv_conc, dv_conc=dv_conc, dna_conc=dna_conc)

    # Full primer: overlap (5') + binding (3')
    full_primer = best_overlap + binding_region

    # Warnings
    if len(full_primer) > 60:
        warnings.append(f"Primer is {len(full_primer)}nt — long primers may have synthesis issues")

    if best_overlap_tm.tm < 48:
        warnings.append(f"Overlap Tm {best_overlap_tm.tm}°C is below 48°C — may reduce assembly efficiency")
    elif best_overlap_tm.tm > 65:
        warnings.append(f"Overlap Tm {best_overlap_tm.tm}°C is above 65°C — consider shorter overlap")

    if binding_tm.tm < 55:
        warnings.append(f"Binding Tm {binding_tm.tm}°C is low — may need to extend binding region")

    if abs(best_overlap_tm.tm - overlap_target_tm) > 5:
        warnings.append(
            f"Overlap Tm {best_overlap_tm.tm}°C deviates from target {overlap_target_tm}°C by "
            f"{abs(best_overlap_tm.tm - overlap_target_tm):.1f}°C"
        )

    return GibsonOverlap(
        overlap_seq=best_overlap,
        overlap_tm=best_overlap_tm,
        binding_seq=binding_region,
        binding_tm=binding_tm,
        full_primer=full_primer,
        overlap_length=len(best_overlap),
        warnings=warnings,
    )


# =============================================================================
# GOLDEN GATE ASSEMBLY
# =============================================================================

# High-fidelity 4bp overhang sets from Potapov et al. 2018 (NEB)
# These are pre-validated sets with >95% correct assembly
# Format: list of 4bp overhangs that are mutually compatible
GOLDEN_GATE_FIDELITY_SETS = {
    # 10-fragment assembly set (NEB recommended)
    '10_fragment': [
        'AACG', 'AATG', 'ATAG', 'ATCA', 'GCAA',
        'GCTG', 'GGTA', 'TACG', 'TAGA', 'TGCC',
        'TTCG',  # 11 overhangs = 10 fragments + vector
    ],
    # 4-fragment set (high fidelity)
    '4_fragment': [
        'AATG', 'GCAA', 'TACG', 'TGCC', 'TTCG',
    ],
}

# All 16 possible 4bp reverse-complement palindromes
PALINDROMIC_4BP = {
    'AATT', 'ATAT', 'AGCT', 'ACGT',
    'TATA', 'TTAA', 'TGCA', 'TCGA',
    'GATC', 'GTAC', 'GGCC', 'GCGC',
    'CATG', 'CTAG', 'CGCG', 'CCGG',
}

# Type IIS restriction enzyme cut sites
TYPE_IIS_SITES = {
    'BsaI':  ('GGTCTC', 1),   # cuts 1nt downstream on top, 4nt on bottom
    'BsmBI': ('CGTCTC', 1),
    'BbsI':  ('GAAGAC', 2),
    'SapI':  ('GCTCTTC', 1),
    'BpiI':  ('GAAGAC', 2),
}


@dataclass
class OverhangCheck:
    """Result of Golden Gate overhang analysis."""
    overhang: str
    is_palindromic: bool
    is_in_fidelity_set: bool
    gc_content: float
    reverse_complement: str
    warnings: List[str] = field(default_factory=list)


@dataclass
class GoldenGateResult:
    """Full Golden Gate assembly analysis."""
    overhangs: List[OverhangCheck]
    all_unique: bool
    enzyme: str
    internal_sites: List[Tuple[str, int]]  # (sequence_name, position)
    warnings: List[str] = field(default_factory=list)


def _reverse_complement(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq.upper()))


def check_overhang(overhang: str, fidelity_set: str = '10_fragment') -> OverhangCheck:
    """
    Check a single 4bp overhang for Golden Gate compatibility.
    """
    oh = overhang.upper().strip()
    warnings = []

    if len(oh) != 4:
        warnings.append(f"Overhang must be 4bp, got {len(oh)}")

    is_palindromic = oh in PALINDROMIC_4BP
    if is_palindromic:
        warnings.append(f"{oh} is palindromic — will self-ligate")

    rc = _reverse_complement(oh)
    if oh == rc:
        warnings.append(f"{oh} is self-complementary")

    gc = sum(1 for b in oh if b in 'GC') / 4
    if gc == 0:
        warnings.append("No GC content — weak ligation")
    elif gc == 1.0:
        warnings.append("All GC — may cause secondary structure issues")

    fidelity_overhangs = GOLDEN_GATE_FIDELITY_SETS.get(fidelity_set, [])
    in_fidelity = oh in fidelity_overhangs

    return OverhangCheck(
        overhang=oh,
        is_palindromic=is_palindromic,
        is_in_fidelity_set=in_fidelity,
        gc_content=gc,
        reverse_complement=rc,
        warnings=warnings,
    )


def check_golden_gate(
    overhangs: List[str],
    sequences: Optional[Dict[str, str]] = None,
    primers: Optional[Dict[str, str]] = None,
    enzyme: str = 'BsaI',
    fidelity_set: str = '10_fragment',
) -> GoldenGateResult:
    """
    Validate a set of Golden Gate overhangs.

    Args:
        overhangs: List of 4bp overhang sequences
        sequences: Optional dict of {name: sequence} to scan for internal enzyme sites
        primers: Optional dict of {name: sequence} — primers are also scanned for
                 internal enzyme sites (a common mistake: accidentally engineering
                 a BsaI site into the binding region or spacer)
        enzyme: Type IIS enzyme name
        fidelity_set: Which pre-validated set to check against

    Note on uniqueness: The collision check is a strict 1:1 string match —
    exact sequence identity or reverse-complement identity. It does NOT compute
    mismatch penalties or edit distance. For custom overhangs outside the NEB
    validated sets, near-misses (1bp difference) may still cause reduced
    fidelity; use the NEB Ligation Fidelity Viewer for edge cases.
    """
    warnings = []

    # Check each overhang
    oh_results = [check_overhang(oh, fidelity_set) for oh in overhangs]

    # Check uniqueness (strict string match + reverse complement)
    oh_seqs = [r.overhang for r in oh_results]
    rc_seqs = [r.reverse_complement for r in oh_results]
    all_seqs = oh_seqs + rc_seqs

    all_unique = len(set(all_seqs)) == len(all_seqs)
    if not all_unique:
        # Find the duplicates
        seen = set()
        for s in all_seqs:
            if s in seen:
                warnings.append(f"Duplicate or RC collision: {s}")
            seen.add(s)

    # Scan for internal enzyme sites in part sequences AND primers
    internal_sites = []
    if enzyme in TYPE_IIS_SITES:
        site_seq, _ = TYPE_IIS_SITES[enzyme]
        site_rc = _reverse_complement(site_seq)

        # Combine sequences and primers for scanning
        all_seqs_to_scan = {}
        if sequences:
            all_seqs_to_scan.update(sequences)
        if primers:
            all_seqs_to_scan.update(
                {f"primer:{k}": v for k, v in primers.items()}
            )

        for name, seq in all_seqs_to_scan.items():
            seq_upper = seq.upper()
            for pos in range(len(seq_upper) - len(site_seq) + 1):
                window = seq_upper[pos:pos + len(site_seq)]
                if window == site_seq or window == site_rc:
                    internal_sites.append((name, pos))
                    if name.startswith('primer:'):
                        warnings.append(
                            f"Internal {enzyme} site in {name} at position {pos} — "
                            f"primer will cleave itself during Golden Gate reaction"
                        )
                    else:
                        warnings.append(f"Internal {enzyme} site in {name} at position {pos}")

    # Overall assessment
    palindromic_count = sum(1 for r in oh_results if r.is_palindromic)
    if palindromic_count:
        warnings.append(f"{palindromic_count} palindromic overhang(s) — expect self-ligation")

    non_fidelity = sum(1 for r in oh_results if not r.is_in_fidelity_set)
    if non_fidelity:
        warnings.append(
            f"{non_fidelity} overhang(s) not in NEB high-fidelity set — "
            f"consider using validated overhangs"
        )

    return GoldenGateResult(
        overhangs=oh_results,
        all_unique=all_unique,
        enzyme=enzyme,
        internal_sites=internal_sites,
        warnings=warnings,
    )


# =============================================================================
# ASSEMBLY METHOD SELECTOR
# =============================================================================

@dataclass
class AssemblyRecommendation:
    """Recommendation for which assembly method to use."""
    method: str          # 'gibson', 'golden_gate', 'restriction_ligation', 'overlap_pcr'
    reason: str
    fragment_count: int
    total_size_bp: int
    warnings: List[str] = field(default_factory=list)


def recommend_assembly(
    fragment_count: int,
    total_size_bp: int,
    needs_scarless: bool = True,
    has_repeats: bool = False,
) -> AssemblyRecommendation:
    """
    Recommend an assembly method based on project parameters.
    """
    warnings = []

    if fragment_count <= 2 and not needs_scarless:
        method = 'restriction_ligation'
        reason = "Simple 2-fragment join — restriction/ligation is cheapest and most reliable"
    elif fragment_count <= 3 and total_size_bp < 10000:
        method = 'gibson'
        reason = "2–3 fragments under 10kb — Gibson is fast, scarless, and well-validated"
    elif fragment_count <= 6:
        method = 'gibson'
        reason = "Multi-fragment Gibson works well up to 6 fragments"
        if fragment_count > 4:
            warnings.append("Consider Golden Gate for >4 fragments — higher efficiency")
    elif fragment_count <= 12:
        method = 'golden_gate'
        reason = "6+ fragments — Golden Gate gives higher efficiency with validated overhangs"
    else:
        method = 'golden_gate'
        reason = "Large assembly — Golden Gate with hierarchical strategy"
        warnings.append("Consider splitting into sub-assemblies for >12 fragments")

    if has_repeats and method == 'gibson':
        warnings.append(
            "Repeated sequences may cause misassembly with Gibson — "
            "consider Golden Gate or overlap extension PCR"
        )

    if total_size_bp > 20000:
        warnings.append("Large construct — verify with restriction digest after assembly")

    return AssemblyRecommendation(
        method=method,
        reason=reason,
        fragment_count=fragment_count,
        total_size_bp=total_size_bp,
        warnings=warnings,
    )
