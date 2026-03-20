"""
primerdesignr.mathews_hairpin — Mathews/Turner RNA parameter hairpin heuristic

A secondary hairpin scanner providing a heuristic second opinion on
ssDNA hairpin prediction. SantaLucia DNA parameters include a terminal
A-T penalty on hairpin closing pairs; Mathews/Turner RNA parameters do
not. Binet et al. (2023) found that RNA parameters predicted ssDNA
secondary structure slightly more accurately (49% vs 43% identical to
experiment) in mfold, partly attributable to this penalty difference.

However, that study evaluated full secondary-structure prediction on
structured aptamers, not targeted hairpin scanning on short primers.
This module is therefore a heuristic flag for manual review, not a
ground-truth engine. It is most useful when SantaLucia and Mathews
disagree on a specific hairpin with an A-T closing pair.

This is NOT a full RNA folding engine. It scans for hairpin structures
with stem ≥3 bp and loop 3-30 nt using Mathews/Turner NN stacking
parameters.

References:
    Mathews DH et al. (1999) J Mol Biol 288:911-940
    Xia T et al. (1998) Biochemistry 37:14719-14735
    Binet et al. (2023) BMC Bioinformatics 24:422
"""

import math
from typing import List, Tuple

# Turner/Mathews RNA nearest-neighbor parameters
# ΔH (kcal/mol), ΔS (cal/mol·K)
# We map T→U for the lookup since these are RNA params applied to DNA
_NN = {
    'AA/UU': (-6.8, -19.0),
    'AU/UA': (-9.4, -26.7),
    'UA/AU': (-7.7, -20.5),
    'CA/GU': (-10.4, -26.9),
    'GU/CA': (-10.2, -26.2),
    'CU/GA': (-7.6, -19.2),
    'GA/CU': (-12.4, -32.5),
    'CG/GC': (-10.6, -26.7),
    'GC/CG': (-14.9, -36.9),
    'GG/CC': (-13.4, -32.7),
}

# Hairpin loop initiation ΔG at 37°C (kcal/mol)
_LOOP_DG = {
    3: 5.6,
    4: 4.3,
    5: 4.6,
    6: 4.4,
    7: 4.6,
    8: 4.8,
    9: 4.9,
    10: 4.9,
}

_COMP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def _to_rna(base: str) -> str:
    return 'U' if base == 'T' else base


def _nn_dg(top: str, bottom: str, temp_k: float) -> float:
    """Look up NN ΔG for a dinucleotide pair, converting DNA→RNA keys."""
    top_r = _to_rna(top[0]) + _to_rna(top[1])
    bot_r = _to_rna(bottom[0]) + _to_rna(bottom[1])

    key = f"{top_r}/{bot_r}"
    if key in _NN:
        dH, dS = _NN[key]
    else:
        # Try reverse complement orientation
        rev_key = f"{bot_r[::-1]}/{top_r[::-1]}"
        if rev_key in _NN:
            dH, dS = _NN[rev_key]
        else:
            dH, dS = -8.0, -22.0  # Fallback average

    return dH - (temp_k * dS / 1000)


def _loop_dg(loop_size: int, temp_k: float) -> float:
    """Loop initiation penalty."""
    if loop_size in _LOOP_DG:
        return _LOOP_DG[loop_size]
    # Jacobson-Stockmayer extrapolation for large loops
    R = 1.987 / 1000  # kcal/mol·K
    return _LOOP_DG[10] + 2.44 * R * temp_k * math.log(loop_size / 10)


def _find_all_hairpins(seq: str, min_stem: int = 3, min_loop: int = 3,
                       max_loop: int = 10) -> List[Tuple[List[Tuple[int, int]], int]]:
    """
    Scan sequence for all possible hairpin configurations.

    Returns list of (stem_pairs, loop_size) tuples.
    stem_pairs: [(i, j), ...] where i < j, ordered from closing pair outward.
    """
    n = len(seq)
    hairpins = []

    for loop_start in range(min_stem, n - min_stem - min_loop + 1):
        for loop_size in range(min_loop, min(max_loop + 1, n - loop_start - min_stem + 1)):
            loop_end = loop_start + loop_size - 1

            # Extend stem outward from loop
            pairs = []
            i = loop_start - 1
            j = loop_end + 1

            while i >= 0 and j < n:
                if _COMP.get(seq[i]) == seq[j]:
                    pairs.append((i, j))
                    i -= 1
                    j += 1
                else:
                    break

            if len(pairs) >= min_stem:
                hairpins.append((pairs, loop_size))

    return hairpins


def calc_hairpin_dg(seq: str, temp: float = 37.0) -> float:
    """
    Calculate the most stable hairpin ΔG using Mathews/Turner RNA parameters.

    Key difference from SantaLucia: no terminal A-T penalty on closing pairs.

    Args:
        seq: DNA sequence (5'→3')
        temp: Temperature in °C

    Returns:
        Most stable (most negative) hairpin ΔG in kcal/mol.
        Returns 0.0 if no hairpin found.
    """
    seq = seq.upper()
    temp_k = temp + 273.15

    hairpins = _find_all_hairpins(seq, min_stem=3, min_loop=3, max_loop=30)

    if not hairpins:
        return 0.0

    best_dg = 0.0

    for pairs, loop_size in hairpins:
        dg = _loop_dg(loop_size, temp_k)

        # Sum nearest-neighbor stacking in the stem
        for k in range(len(pairs) - 1):
            i1, j1 = pairs[k]      # Closing pair (inner)
            i2, j2 = pairs[k + 1]  # Next pair (outer)

            # 5'→3' on left arm (top): outer→inner = i2, i1
            top = seq[i2] + seq[i1]
            # 3'→5' on right arm (bottom): outer→inner = j2, j1
            bottom = seq[j2] + seq[j1]

            dg += _nn_dg(top, bottom, temp_k)

        # No terminal A-T penalty in Mathews/Turner params.
        # This is the parameter difference Binet et al. identified as
        # potentially affecting ssDNA hairpin predictions at A-T closing pairs.

        if dg < best_dg:
            best_dg = dg

    return round(best_dg, 2)
