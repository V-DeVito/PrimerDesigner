"""
Primer design workflows built on Primer3.

This module keeps design orchestration separate from thermodynamic analysis:
Primer3 proposes primer pairs, then primerdesignr re-analyzes each candidate
with the same warning model used by the analysis workflow.
"""

from dataclasses import dataclass, field
from typing import Optional
import re

import primer3

from .thermo import DimerResult, PairReport, analyze_pair


_COMP = str.maketrans("ATGC", "TACG")


def _reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]


@dataclass
class PrimerCoordinates:
    """Template coordinates for a designed primer."""

    start: int
    end: int
    length: int
    strand: str


@dataclass
class PrimerPairCandidate:
    """A ranked Primer3 candidate plus primerdesignr analysis."""

    rank: int
    pair: PairReport
    forward_coords: PrimerCoordinates
    reverse_coords: PrimerCoordinates
    product_size: int
    primer3_pair_penalty: float
    primer3_pair_compl_any: Optional[float] = None
    primer3_pair_compl_end: Optional[float] = None
    explanations: list[str] = field(default_factory=list)

    @property
    def pair_report(self) -> PairReport:
        return self.pair

    @property
    def heterodimer(self) -> DimerResult:
        return self.pair.heterodimer

    @property
    def warnings(self) -> list[str]:
        return self.pair.warnings + self.pair.forward.warnings + self.pair.reverse.warnings


@dataclass
class PrimerDesignResult:
    """Result of a PCR primer design request."""

    template_length: int
    product_size_range: tuple[int, int]
    target_region: Optional[tuple[int, int]]
    candidates: list[PrimerPairCandidate]
    primer3_explain: dict[str, str]
    warnings: list[str] = field(default_factory=list)


def _clean_template(template: str) -> str:
    lines = []
    for line in template.splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith(">"):
            lines.append(stripped)
    seq = re.sub(r"\s+", "", "".join(lines)).upper()
    if re.search(r"[^ATGC]", seq):
        raise ValueError("Template contains non-ATGC characters")
    return seq


def _coord_from_left(value: list[int]) -> PrimerCoordinates:
    start, length = value
    return PrimerCoordinates(start=start, end=start + length - 1, length=length, strand="+")


def _coord_from_right(value: list[int]) -> PrimerCoordinates:
    end, length = value
    return PrimerCoordinates(start=end - length + 1, end=end, length=length, strand="-")


def _candidate_explanations(
    pair: PairReport,
    product_size: int,
    penalty: float,
    design_mode: str,
) -> list[str]:
    explanations = [
        f"Product size is {product_size} bp.",
        f"Forward/reverse Tm difference is {pair.tm_difference:.1f} C.",
        f"Pair heterodimer delta G is {pair.heterodimer.dg:.2f} kcal/mol.",
    ]
    if design_mode == "exact":
        explanations.append("Exact-sequence mode fixes primer binding to the submitted 5' and 3' ends, then varies primer length.")
    else:
        explanations.append(f"Primer3 pair penalty is {penalty:.2f}; lower values are closer to the configured amplicon targets.")

    if pair.tm_difference <= 2:
        explanations.append("Tm values are tightly matched.")
    elif pair.tm_difference <= 5:
        explanations.append("Tm values are acceptable but not tightly matched.")
    else:
        explanations.append("Tm mismatch exceeds the default 5 C warning threshold.")

    if not pair.warnings and not pair.forward.warnings and not pair.reverse.warnings:
        explanations.append("No PrimerDesigner warnings were raised for this pair.")

    return explanations


def _dimer_excess(candidate: PrimerPairCandidate) -> float:
    """Amount by which dimers exceed PrimerDesigner's -9 kcal/mol warning line."""
    values = [
        candidate.heterodimer.dg,
        candidate.pair.forward.homodimer.dg,
        candidate.pair.reverse.homodimer.dg,
    ]
    return round(sum(max(0.0, -dg - 9.0) for dg in values), 3)


def _hairpin_excess(candidate: PrimerPairCandidate) -> float:
    values = [
        candidate.pair.forward.hairpin.worst_dg,
        candidate.pair.reverse.hairpin.worst_dg,
    ]
    return round(sum(max(0.0, -dg - 3.0) for dg in values), 3)


def _continuous_quality(candidate: PrimerPairCandidate) -> float:
    pair = candidate.pair
    primers = [pair.forward, pair.reverse]
    tm_score = sum(abs(primer.tm.tm - 60.0) / 6.0 for primer in primers)
    gc_score = sum(abs((primer.gc_content * 100) - 50.0) / 15.0 for primer in primers)
    length_score = sum(abs(primer.length - 20) / 8.0 for primer in primers)
    hairpin_score = sum(max(0.0, -primer.hairpin.worst_dg) / 8.0 for primer in primers)
    dimer_score = (
        max(0.0, -pair.heterodimer.dg) +
        max(0.0, -pair.forward.homodimer.dg) +
        max(0.0, -pair.reverse.homodimer.dg)
    ) / 18.0

    return round(
        (pair.tm_difference * 0.9) +
        (tm_score * 0.7) +
        (gc_score * 0.35) +
        (length_score * 0.25) +
        (hairpin_score * 0.45) +
        (dimer_score * 0.75),
        3,
    )


def _danger_excess(candidate: PrimerPairCandidate) -> float:
    return round(_dimer_excess(candidate) + _hairpin_excess(candidate), 3)


def _candidate_score(candidate: PrimerPairCandidate) -> tuple[float, float, float, float]:
    """
    Rank by PrimerDesigner's launch-facing quality first, then Primer3 penalty.

    Primer3 is excellent at finding legal primer pairs, but its pair penalty
    does not map one-to-one to the warnings this app explains to users. We
    therefore request a deeper Primer3 pool and promote candidates with fewer
    PrimerDesigner warnings.
    """
    return (
        len(candidate.warnings),
        _danger_excess(candidate),
        _continuous_quality(candidate),
        candidate.primer3_pair_penalty,
    )


def _near_duplicate(a: PrimerPairCandidate, b: PrimerPairCandidate) -> bool:
    """True when two pairs represent the same practical primer choice."""
    return (
        abs(a.forward_coords.start - b.forward_coords.start) <= 8
        and abs(a.forward_coords.end - b.forward_coords.end) <= 8
        and abs(a.reverse_coords.start - b.reverse_coords.start) <= 8
        and abs(a.reverse_coords.end - b.reverse_coords.end) <= 8
        and abs(a.product_size - b.product_size) <= 16
    )


def _select_diverse_candidates(
    candidates: list[PrimerPairCandidate],
    primer_count: int,
) -> tuple[list[PrimerPairCandidate], int]:
    """
    Keep the best pair from each local coordinate cluster.

    Primer3 often returns many one- or two-base shifts around the same binding
    sites. Those are useful internally for scoring, but they are not meaningful
    alternatives for a user-facing candidate list.
    """
    selected: list[PrimerPairCandidate] = []
    hidden_near_duplicates = 0

    for candidate in sorted(candidates, key=_candidate_score):
        if any(_near_duplicate(candidate, kept) for kept in selected):
            hidden_near_duplicates += 1
            continue
        selected.append(candidate)
        if len(selected) >= primer_count:
            break

    return selected, hidden_near_duplicates


def _exact_end_candidates(
    seq: str,
    primer_count: int,
    mv_conc: float,
    dv_conc: float,
    dntp_conc: float,
    dna_conc: float,
) -> list[PrimerPairCandidate]:
    candidates = []
    for forward_len in range(17, min(30, len(seq) - 17) + 1):
        fwd_seq = seq[:forward_len]
        for reverse_len in range(17, min(30, len(seq) - forward_len) + 1):
            rev_seq = _reverse_complement(seq[-reverse_len:])
            pair = analyze_pair(
                fwd_seq,
                rev_seq,
                mv_conc=mv_conc,
                dv_conc=dv_conc,
                dntp_conc=dntp_conc,
                dna_conc=dna_conc,
            )
            candidate = PrimerPairCandidate(
                rank=0,
                pair=pair,
                forward_coords=PrimerCoordinates(start=0, end=forward_len - 1, length=forward_len, strand="+"),
                reverse_coords=PrimerCoordinates(
                    start=len(seq) - reverse_len,
                    end=len(seq) - 1,
                    length=reverse_len,
                    strand="-",
                ),
                product_size=len(seq),
                primer3_pair_penalty=0.0,
            )
            candidate.explanations = _candidate_explanations(pair, len(seq), 0.0, "exact")
            candidates.append(candidate)

    return sorted(candidates, key=_candidate_score)[:primer_count]


def design_pcr_primers(
    template: str,
    product_min: int = 120,
    product_max: int = 500,
    primer_count: int = 5,
    target_start: Optional[int] = None,
    target_length: Optional[int] = None,
    mv_conc: float = 50.0,
    dv_conc: float = 0.0,
    dntp_conc: float = 0.0,
    dna_conc: float = 250.0,
    primer_min_tm: float = 57.0,
    primer_opt_tm: float = 60.0,
    primer_max_tm: float = 63.0,
    design_mode: str = "exact",
) -> PrimerDesignResult:
    """
    Design PCR primer pairs for a DNA template.

    Coordinates are 0-based and inclusive. If target_start/target_length are
    supplied, Primer3 is asked to include that region in the amplicon.
    """
    seq = _clean_template(template)
    target_region = None
    sequence_args: dict[str, object] = {
        "SEQUENCE_ID": "template",
        "SEQUENCE_TEMPLATE": seq,
    }

    if target_start is not None and target_length is not None:
        target_region = (target_start, target_length)
        sequence_args["SEQUENCE_TARGET"] = [target_start, target_length]

    if design_mode == "exact":
        candidates = _exact_end_candidates(
            seq,
            1,
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
        )
        warnings = []
        if candidates and all(candidate.warnings for candidate in candidates):
            warnings.append(
                "No warning-free exact-end primer pair was found. Provide upstream/downstream "
                "flanking sequence or allow an internal amplicon so PrimerDesigner has more "
                "binding sites to consider."
            )

        for rank, candidate in enumerate(candidates, start=1):
            candidate.rank = rank

        return PrimerDesignResult(
            template_length=len(seq),
            product_size_range=(len(seq), len(seq)),
            target_region=target_region,
            candidates=candidates,
            primer3_explain={},
            warnings=warnings,
        )

    internal_count = min(100, max(primer_count, primer_count * 12, 30))
    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_NUM_RETURN": internal_count,
        "PRIMER_PRODUCT_SIZE_RANGE": [[product_min, product_max]],
        "PRIMER_MIN_SIZE": 17,
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MAX_SIZE": 30,
        "PRIMER_MIN_TM": min(primer_min_tm, 55.0),
        "PRIMER_OPT_TM": primer_opt_tm,
        "PRIMER_MAX_TM": max(primer_max_tm, 65.0),
        "PRIMER_MIN_GC": 30,
        "PRIMER_MAX_GC": 70,
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_MAX_SELF_ANY_TH": 45,
        "PRIMER_MAX_SELF_END_TH": 35,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45,
        "PRIMER_PAIR_MAX_COMPL_END_TH": 35,
        "PRIMER_SALT_MONOVALENT": mv_conc,
        "PRIMER_SALT_DIVALENT": dv_conc,
        "PRIMER_DNTP_CONC": dntp_conc,
        "PRIMER_DNA_CONC": dna_conc,
    }

    raw = primer3.bindings.design_primers(sequence_args, global_args)
    returned = int(raw.get("PRIMER_PAIR_NUM_RETURNED", 0))

    candidates: list[PrimerPairCandidate] = []
    for i in range(returned):
        fwd_seq = raw[f"PRIMER_LEFT_{i}_SEQUENCE"]
        rev_seq = raw[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        pair = analyze_pair(
            fwd_seq,
            rev_seq,
            mv_conc=mv_conc,
            dv_conc=dv_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
        )
        product_size = int(raw[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"])
        penalty = float(raw.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0))
        candidate = PrimerPairCandidate(
            rank=i + 1,
            pair=pair,
            forward_coords=_coord_from_left(raw[f"PRIMER_LEFT_{i}"]),
            reverse_coords=_coord_from_right(raw[f"PRIMER_RIGHT_{i}"]),
            product_size=product_size,
            primer3_pair_penalty=round(penalty, 3),
            primer3_pair_compl_any=raw.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
            primer3_pair_compl_end=raw.get(f"PRIMER_PAIR_{i}_COMPL_END_TH"),
        )
        candidate.explanations = _candidate_explanations(pair, product_size, penalty, "amplicon")
        candidates.append(candidate)

    explain = {
        key: str(value)
        for key, value in raw.items()
        if key.endswith("_EXPLAIN") or key in {"PRIMER_WARNING", "PRIMER_ERROR"}
    }

    warnings = []
    if not candidates:
        warnings.append(
            "Primer3 did not return candidates for the selected constraints. "
            "Relax product size or target constraints, or provide more upstream/downstream sequence "
            "so PrimerDesigner has more binding sites to consider."
        )

    candidates, hidden_near_duplicates = _select_diverse_candidates(candidates, primer_count)
    for rank, candidate in enumerate(candidates, start=1):
        candidate.rank = rank

    if hidden_near_duplicates and len(candidates) < primer_count:
        warnings.append(
            f"{hidden_near_duplicates} near-duplicate primer pair(s) were hidden; fewer "
            "meaningfully distinct candidates were available under the current constraints."
        )

    if candidates and all(candidate.warnings for candidate in candidates):
        warnings.append(
            "No warning-free primer pair was found under the current template and constraints. "
            "Relax product size or target constraints, or provide more upstream/downstream sequence "
            "so PrimerDesigner has more binding sites to consider."
        )

    return PrimerDesignResult(
        template_length=len(seq),
        product_size_range=(product_min, product_max),
        target_region=target_region,
        candidates=candidates,
        primer3_explain=explain,
        warnings=warnings,
    )
