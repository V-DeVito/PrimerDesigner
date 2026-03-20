# primerdesignr — Technical Methods Guide

primerdesignr is a pragmatic primer-screening tool that combines Primer3 thermodynamic analysis with a secondary ssDNA hairpin heuristic to catch edge cases that single-engine workflows may miss.

---

## 1. Melting Temperature (Tm)

**Engine:** primer3-py (`calc_tm`)

**What it is:** The temperature at which 50% of primer-template duplexes are dissociated.

**Method:** The nearest-neighbor (NN) model treats duplex stability as the sum of stacking interactions between adjacent base pairs. Each dinucleotide step contributes an experimentally measured enthalpy (ΔH) and entropy (ΔS). The total values are:

    ΔH_total = ΔH_initiation + Σ ΔH_i
    ΔS_total = ΔS_initiation + Σ ΔS_i + R·ln(Ct/4)
    Tm = ΔH_total / ΔS_total

Where Ct is the total strand concentration and R = 1.987 cal/mol·K. Initiation parameters account for the nucleation cost at each duplex terminus and differ for G-C vs A-T terminal pairs. The specific initiation constants used are those implemented in primer3-py's SantaLucia engine; consult the primer3-py source or SantaLucia & Hicks (2004) Table 2 for exact values.

**Salt correction:** SantaLucia (1998) monovalent correction by default. When Mg2+ is present and dominates, the Owczarzy et al. (2008) correction applies.

**primerdesignr default conditions vs primer3-py defaults:**

primerdesignr overrides primer3-py's library defaults to match IDT OligoAnalyzer conventions:

| Parameter | primer3-py default | primerdesignr default | Notes |
|---|---|---|---|
| Na+ (mv_conc) | 50.0 mM | 50.0 mM | Same |
| Mg2+ (dv_conc) | 1.5 mM | 0.0 mM | **Override** — IDT default is 0 |
| dNTP (dntp_conc) | 0.6 mM | 0.0 mM | **Override** — IDT default is 0 |
| Primer (dna_conc) | 50.0 nM | 250.0 nM | **Override** — IDT default is 250 |
| Salt correction | santalucia | santalucia | Same |
| Tm method | santalucia | santalucia | Same |

These overrides are applied consistently across all calculations: Tm, homodimer, heterodimer, and 3' end stability. For realistic PCR conditions (e.g., 1× Taq buffer), pass `dv_conc=1.5, dntp_conc=0.2`. This typically increases Tm by 6–7°C. All functions accept user-specified conditions.

**Parameter source:** SantaLucia Jr, J. & Hicks, D. (2004) Annu Rev Biophys Biomol Struct 33:415–440.

---

## 2. Hairpin ΔG (Dual Engine)

**Engines:** seqfold (SantaLucia 2004 DNA, Zuker DP) + custom module (Mathews/Turner RNA parameters)

**What it is:** The minimum free energy (MFE) of the most stable intramolecular hairpin structure the primer can form by folding back on itself. More negative ΔG = more stable = more problematic.

**Why two engines:**

SantaLucia DNA parameters include an energetic penalty for hairpin stems that close with A-T base pairs. The Mathews/Turner RNA parameters do not. Binet et al. (BMC Bioinformatics, 2023) evaluated full secondary-structure prediction tools on a dataset of experimentally determined ssDNA structures and found that mfold with RNA (default) parameters predicted 49% of ssDNA structures identically to experiment vs 43% with the SantaLucia DNA model. The authors specifically noted that the terminal A-T penalties in the DNA model may contribute to this difference in some cases.

However, that study evaluated complete secondary-structure prediction on structured aptamers (40–200+ nt), not targeted hairpin scanning on short PCR primers. The same paper also concludes more broadly that DNA thermodynamic models did not generally improve ssDNA structure prediction overall. Our Mathews engine should therefore be understood as a **secondary heuristic view** for ssDNA hairpin prediction — not as a more accurate ground truth. It is most useful as a flag for manual review when SantaLucia and Mathews disagree on a specific hairpin, particularly one with an A-T closing pair.

**Engine 1 — seqfold (SantaLucia 2004):**

seqfold implements the Zuker (1981) dynamic programming algorithm for MFE structure prediction with energy functions from SantaLucia & Hicks (2004). It produces dot-bracket notation and accounts for base-pair stacking, hairpin loop initiation, internal loops, bulges, and multiloop penalties.

**Engine 2 — mathews_hairpin (Turner/Mathews RNA parameters):**

A targeted hairpin scanner using Turner/Mathews (1999) nearest-neighbor stacking parameters with T→U mapping. It scans for hairpin structures with stem lengths ≥3 bp and loop sizes 3–30 nt. The key difference from Engine 1: no terminal A-T penalty on closing pairs. This is the specific parameter difference that Binet et al. identified as potentially affecting some ssDNA hairpin predictions.

**Threshold:** ΔG < -3.0 kcal/mol flags a hairpin concern. This is a tool-level screening threshold under the stated thermodynamic conditions, not an absolute biochemical boundary — the actual stability depends on salt concentration, temperature, and primer context.

**Disagreement warning:** When SantaLucia reports safe (ΔG > -3.0) but Mathews reports problematic (ΔG < -3.0), the warning reads: "May be under-penalized by DNA parameters when stem closes with A-T; review manually." This framing reflects the literature without overstating the Mathews result.

---

## 3. Homodimer ΔG (Self-Dimer)

**Engine:** primer3-py (`calc_homodimer`)

**What it is:** The ΔG of the most stable structure formed when two copies of the same primer hybridize to each other in an antiparallel orientation.

**Method:** Primer3's thermodynamic alignment engine evaluates all possible alignments of a primer against its reverse complement. The engine uses nearest-neighbor thermodynamic parameters to score alignments, accounting for Watson-Crick base pairs, mismatches, gaps, bulges, and loop structures — not just contiguous complementary stretches. It returns the most stable (most negative) ΔG found across all alignments.

**Threshold:** ΔG < -9.0 kcal/mol is flagged. This is a primerdesignr screening threshold under the tool's default salt/temperature conditions, consistent with common industry practice. It is not an absolute biochemical boundary — changing salt concentration or temperature will shift the actual ΔG.

**3' end stability:** In addition to overall ΔG, primerdesignr evaluates the 3' end alignment using primer3-py's `calc_end_stability()`. A dimer with a modest overall ΔG can ruin a PCR if the complementarity is at the 3' end, because the polymerase extends from the 3' terminus. These "weak" but well-positioned dimers act as primers for each other, generating primer-dimer artifacts. The 3' end stability threshold is -5.0 kcal/mol.

**Conditions:** primerdesignr passes its own default conditions (see Section 1 table) to all primer3-py thermo calls, overriding primer3-py's library defaults. User-specified conditions propagate consistently through all calculations.

---

## 4. Heterodimer ΔG (Cross-Dimer)

**Engine:** primer3-py (`calc_heterodimer`)

**What it is:** The ΔG of the most stable structure formed when two different primers hybridize to each other. Same thermodynamic alignment method as homodimer analysis (Section 3), evaluated across all alignments of primer 1 against the reverse complement of primer 2. Same -9.0 kcal/mol screening threshold.

**3' end stability:** Both orientations are checked (seq1's 3' on seq2, and seq2's 3' on seq1) via `calc_end_stability()`. The worse (more negative) value is reported. Same -5.0 kcal/mol threshold.

**Cross-dimer matrix:** For a set of N primers, primerdesignr computes all N(N-1)/2 unique pairwise combinations.

---

## 5. GC Content and 3' GC Clamp

**GC content:** The fraction of G+C bases in the primer. Ideal range: 40–60%.

**3' GC clamp:** primerdesignr checks the last 5 bases for GC content. The ideal GC clamp is 1–2 G/C bases in the final 5, ideally with a G or C at the 3' terminus. This anchors polymerase initiation without making the 3' end so sticky that it tolerates mismatches and misprimes. Graded warnings: 0 G/C (weak anchoring), 3 G/C (acceptable but not ideal), ≥4 G/C (mispriming risk).

---

## 6. Golden Gate Overhang Validation

**What it is:** Golden Gate assembly uses Type IIS restriction enzymes (BsaI, BsmBI, etc.) that cut outside their recognition sequence, leaving user-defined 4-bp overhangs.

**Checks performed:**

**Palindrome detection:** All 16 possible 4-bp reverse-complement palindromes (AATT, ATAT, AGCT, ACGT, TATA, TTAA, TGCA, TCGA, GATC, GTAC, GGCC, GCGC, CATG, CTAG, CGCG, CCGG) are flagged because they are self-complementary and will self-ligate.

**Uniqueness (strict string match):** Each overhang must be distinct from every other overhang and from the reverse complement of every other overhang. This is an exact 1:1 string match — no mismatch penalty or edit distance calculation is performed. For custom overhangs outside the NEB validated sets, near-misses (1bp difference) may still cause reduced fidelity; consult the NEB Ligation Fidelity Viewer for edge cases.

**NEB high-fidelity sets:** Potapov et al. (ACS Synth Biol, 2018) experimentally determined overhang combinations with >95% correct assembly. primerdesignr checks against the NEB-validated 10-fragment and 4-fragment sets.

**Internal site scanning (parts and primers):** Both part sequences and primer sequences are scanned for internal occurrences of the Type IIS recognition site (and its reverse complement). Scanning primers is important because it is common to accidentally engineer a restriction site into a primer's binding region or spacer, causing the amplified product to cleave itself during the Golden Gate reaction.

---

## 7. Gibson Overlap Tm

**What it is:** In Gibson assembly, adjacent fragments share a 15–40 bp overlap. The overlap Tm is calculated using the same SantaLucia nearest-neighbor method as primer Tm (Section 1). primerdesignr scans overlap lengths from 15 to 40 bp and selects the length whose Tm is closest to the target (default 50°C per NEB guidelines).

---

## Warnings

All ΔG thresholds listed below are **primerdesignr screening thresholds under the tool's default thermodynamic conditions** (50 mM Na+, 0 mM Mg2+, 0 mM dNTP, 250 nM primer, 37°C). They are pragmatic cutoffs consistent with industry practice, not absolute biochemical boundaries. Changing salt concentration, temperature, or primer concentration will shift the actual ΔG values.

| Warning | Condition | Rationale |
|---|---|---|
| Short primer | <18 nt | Too short for specific binding in most genomes |
| Long primer | >30 nt | May reduce specificity; synthesis cost increases |
| Low GC | <40% | Weak binding; low Tm |
| High GC | >60% | Secondary structure risk; mispriming |
| Low Tm | <52°C | May not anneal at standard PCR temperatures |
| High Tm | >68°C | Consider shortening primer |
| Stable hairpin | ΔG < -3.0 kcal/mol (either engine) | Primer self-folds, blocking template binding |
| Hairpin disagreement | SL safe, MW problematic | May be under-penalized by DNA parameters at A-T closing pair; review manually |
| Stable self-dimer | ΔG < -9.0 kcal/mol | Primers sequester each other |
| 3' self-dimer risk | 3' end ΔG < -5.0 kcal/mol | Polymerase extends from 3', creating primer-dimer artifacts |
| Stable cross-dimer | ΔG < -9.0 kcal/mol | Forward/reverse bind each other |
| 3' cross-dimer risk | 3' end ΔG < -5.0 kcal/mol | Primers can extend off each other at 3' ends |
| No 3' GC | 0 G/C in last 5 bases | Weak 3' anchoring |
| Borderline 3' GC | 3 G/C in last 5 bases | Acceptable but prefer 1–2 for ideal clamp |
| Heavy 3' GC | ≥4 G/C in last 5 bases | Mispriming risk from overly sticky 3' end |
| Homopolymer run | ≥4 consecutive identical bases | Polymerase slippage risk |
| Tm difference | >5°C between pair members | Unequal amplification efficiency |
| Internal Type IIS (primer) | BsaI/BsmBI/etc. in primer | Primer will cleave itself during Golden Gate |
| Internal Type IIS (part) | BsaI/BsmBI/etc. in part | Part cut at unintended position |

---

## Parameter Sources

| Parameter | Source |
|---|---|
| DNA NN stacking ΔH, ΔS | SantaLucia & Hicks (2004) Annu Rev Biophys 33:415 |
| Monovalent salt correction | SantaLucia (1998) PNAS 95:1460 |
| Mg2+ salt correction | Owczarzy et al. (2008) Biochemistry 47:5336 |
| RNA NN stacking parameters | Mathews et al. (1999) J Mol Biol 288:911; Xia et al. (1998) Biochemistry 37:14719 |
| Zuker DP algorithm | Zuker & Stiegler (1981) Nucleic Acids Res 9:133 |
| SL vs Mathews accuracy | Binet et al. (2023) BMC Bioinformatics 24:422 |
| Golden Gate fidelity sets | Potapov et al. (2018) ACS Synth Biol 7:2665 |
| Gibson overlap guidelines | Gibson et al. (2009) Nature Methods 6:343 |

---

## Software Dependencies

| Library | Version | License | Role |
|---|---|---|---|
| primer3-py | ≥2.0 | GPL v2 | Tm, homodimer, heterodimer, 3' end stability |
| seqfold | ≥0.7 | MIT | Hairpin MFE + dot-bracket (SantaLucia 2004) |
| mathews_hairpin (built-in) | — | MIT | Hairpin heuristic second opinion (Turner/Mathews) |

**Licensing note:** primer3-py is GPL v2. primerdesignr uses it as a server-side dependency behind an HTTP API. If distributing primerdesignr as a bundled library or commercial product, the GPL implications for primer3-py should be evaluated by counsel.
