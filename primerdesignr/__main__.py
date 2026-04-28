"""
primerdesignr CLI — Primer analysis from the command line.

Usage:
    python -m primerdesignr analyze SEQ [SEQ...]
    python -m primerdesignr batch primers.txt
    python -m primerdesignr pair FORWARD REVERSE
    python -m primerdesignr golden-gate OVERHANG [OVERHANG...]

Examples:
    python -m primerdesignr analyze GTCTTCACATCGGTTTGAAAGGAGG
    python -m primerdesignr batch my_primers.txt
    python -m primerdesignr pair GTCTTCACATCGGTTTGAAAGGAGG AACCCGCTCCGATTAAAGCTACTTT
    python -m primerdesignr golden-gate AACG AATG ATAG GCAA
"""

import argparse
import sys
import csv
import io
from typing import Dict

from .thermo import (
    analyze_primer, analyze_pair, cross_dimer_matrix,
    PrimerReport, PairReport, DimerResult,
)
from .assembly import check_golden_gate, check_overhang
from .design import design_pcr_primers


def _parse_primer_input(text: str) -> Dict[str, str]:
    """
    Parse primer input in flexible formats:
        NAME SEQUENCE
        NAME\tSEQUENCE
        NAME,SEQUENCE
        >NAME\\nSEQUENCE  (FASTA)
        SEQUENCE  (auto-named)
    """
    primers = {}
    lines = text.strip().split('\n')
    counter = 1

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        # FASTA format
        if line.startswith('>'):
            name = line[1:].strip()
            i += 1
            if i < len(lines):
                seq = lines[i].strip().upper()
                primers[name] = seq
            i += 1
            continue

        # Tab/space/comma separated
        for sep in ['\t', ',', ' ']:
            if sep in line:
                parts = line.split(sep, 1)
                name = parts[0].strip()
                seq = parts[1].strip().upper()
                # Verify it looks like DNA
                if all(c in 'ATGCNRYSWKM' for c in seq):
                    primers[name] = seq
                    break
        else:
            # Just a sequence
            seq = line.upper()
            if all(c in 'ATGCNRYSWKM' for c in seq) and len(seq) >= 10:
                primers[f'Primer_{counter}'] = seq
                counter += 1

        i += 1

    return primers


def _print_report(name: str, r: PrimerReport):
    """Print a single primer report."""
    status = '✓' if not r.warnings else '⚠'
    print(f"\n  {status} {name}: {r.sequence}")
    print(f"    Length: {r.length}bp | GC: {r.gc_content:.0%} | Tm: {r.tm.tm}°C")
    print(f"    Hairpin:    SL={r.hairpin.dg_santalucia:>6.2f}  MW={r.hairpin.dg_mathews:>6.2f} kcal/mol")
    print(f"    Self-dimer: {r.homodimer.dg:>6.2f} kcal/mol")

    if r.hairpin.dot_bracket and any(c in r.hairpin.dot_bracket for c in '()'):
        print(f"    Structure:  {r.hairpin.dot_bracket}")

    if r.hairpin.engines_disagree:
        print(f"    ⚠ SantaLucia/Mathews disagree — likely AT-closing hairpin stem")

    for w in r.warnings:
        print(f"    ⚠ {w}")


def _print_matrix(primers: Dict[str, str], matrix: Dict, names: list):
    """Print cross-dimer matrix."""
    # Column width
    cw = max(len(n) for n in names) + 1
    cw = max(cw, 8)

    print(f"\n  {'':>{cw}}", end='')
    for n in names:
        print(f"{n:>{cw}}", end='')
    print()

    for i, n1 in enumerate(names):
        print(f"  {n1:>{cw}}", end='')
        for j, n2 in enumerate(names):
            if i == j:
                print(f"{'—':>{cw}}", end='')
            elif i < j:
                key = (n1, n2)
                d = matrix[key]
                flag = '⚠' if d.is_problematic else ' '
                print(f"{d.dg:>{cw-1}.1f}{flag}", end='')
            else:
                key = (n2, n1)
                d = matrix[key]
                flag = '⚠' if d.is_problematic else ' '
                print(f"{d.dg:>{cw-1}.1f}{flag}", end='')
        print()


def _export_csv(primers: Dict[str, str], reports: Dict[str, PrimerReport],
                matrix: Dict, names: list) -> str:
    """Export results as CSV string."""
    output = io.StringIO()
    writer = csv.writer(output)

    # Primer analysis
    writer.writerow([
        'Name', 'Sequence', 'Length', 'GC%', 'Tm (°C)',
        'Hairpin ΔG (SL)', 'Hairpin ΔG (MW)', 'Self-dimer ΔG',
        'Dot-bracket', 'Warnings'
    ])
    for name in names:
        r = reports[name]
        writer.writerow([
            name, r.sequence, r.length, f"{r.gc_content:.1%}", r.tm.tm,
            r.hairpin.dg_santalucia, r.hairpin.dg_mathews, r.homodimer.dg,
            r.hairpin.dot_bracket, '; '.join(r.warnings)
        ])

    writer.writerow([])
    writer.writerow(['Cross-dimer matrix (ΔG kcal/mol)'])

    # Matrix header
    writer.writerow([''] + names)
    for i, n1 in enumerate(names):
        row = [n1]
        for j, n2 in enumerate(names):
            if i == j:
                row.append('—')
            elif i < j:
                row.append(f"{matrix[(n1, n2)].dg:.2f}")
            else:
                row.append(f"{matrix[(n2, n1)].dg:.2f}")
        writer.writerow(row)

    return output.getvalue()


def cmd_analyze(args):
    """Analyze one or more primer sequences."""
    # Collect all sequences
    if args.sequences:
        raw = '\n'.join(args.sequences)
    else:
        print("Paste primer sequences (Ctrl+D or Ctrl+Z to finish):")
        raw = sys.stdin.read()

    primers = _parse_primer_input(raw)
    if not primers:
        print("No valid primer sequences found.")
        sys.exit(1)

    names = list(primers.keys())
    reports = {}

    print(f"\n{'=' * 60}")
    print(f"  PRIMER ANALYSIS — {len(primers)} sequence(s)")
    print(f"  Conditions: {args.na}mM Na+, {args.mg}mM Mg2+, {args.dna}nM primer")
    print(f"{'=' * 60}")

    for name, seq in primers.items():
        r = analyze_primer(seq, mv_conc=args.na, dv_conc=args.mg, dna_conc=args.dna)
        reports[name] = r
        _print_report(name, r)

    # Cross-dimer matrix if multiple primers
    if len(primers) > 1:
        matrix = cross_dimer_matrix(primers)

        print(f"\n{'─' * 60}")
        print(f"  CROSS-DIMER MATRIX")
        _print_matrix(primers, matrix, names)

        # Summary
        problematic = [(k, v) for k, v in matrix.items() if v.is_problematic]
        if problematic:
            print(f"\n  ⚠ {len(problematic)} problematic cross-dimer(s):")
            for (n1, n2), v in problematic:
                print(f"    {n1} × {n2}: {v.dg:.2f} kcal/mol")
    else:
        matrix = {}

    # Tm spread
    tms = [r.tm.tm for r in reports.values()]
    if len(tms) > 1:
        print(f"\n  Tm range: {min(tms):.1f}–{max(tms):.1f}°C (spread: {max(tms)-min(tms):.1f}°C)")

    # CSV export
    if args.csv:
        csv_text = _export_csv(primers, reports, matrix, names)
        with open(args.csv, 'w') as f:
            f.write(csv_text)
        print(f"\n  Exported to {args.csv}")


def cmd_pair(args):
    """Analyze a specific primer pair."""
    result = analyze_pair(
        args.forward, args.reverse,
        mv_conc=args.na, dv_conc=args.mg, dna_conc=args.dna
    )

    print(f"\n{'=' * 60}")
    print(f"  PAIR ANALYSIS")
    print(f"{'=' * 60}")

    _print_report('Forward', result.forward)
    _print_report('Reverse', result.reverse)

    status = '✓' if not result.heterodimer.is_problematic else '⚠'
    print(f"\n  Cross-dimer: {result.heterodimer.dg:.2f} kcal/mol  {status}")
    print(f"  Tm difference: {result.tm_difference}°C")

    for w in result.warnings:
        print(f"  ⚠ {w}")


def cmd_batch(args):
    """Analyze primers from a file."""
    with open(args.file) as f:
        raw = f.read()

    # Inject into analyze
    args.sequences = raw.strip().split('\n')
    cmd_analyze(args)


def cmd_design(args):
    """Design PCR primer pairs from a template sequence."""
    if args.file:
        with open(args.file) as f:
            raw = f.read()
    elif args.template:
        raw = args.template
    else:
        print("Paste template sequence or FASTA (Ctrl+D or Ctrl+Z to finish):")
        raw = sys.stdin.read()

    result = design_pcr_primers(
        raw,
        product_min=args.product_min,
        product_max=args.product_max,
        primer_count=args.count,
        mv_conc=args.na,
        dv_conc=args.mg,
        dna_conc=args.dna,
    )

    print(f"\n{'=' * 60}")
    print(f"  PCR PRIMER DESIGN — {result.template_length}bp template")
    print(f"  Product range: {args.product_min}–{args.product_max}bp")
    print(f"{'=' * 60}")

    if not result.candidates:
        print("  No candidates returned.")
        for key, value in result.primer3_explain.items():
            print(f"  {key}: {value}")
        return

    for c in result.candidates:
        print(f"\n  #{c.rank}  product={c.product_size}bp  penalty={c.primer3_pair_penalty:.2f}")
        print(f"    FWD {c.pair.forward.sequence}  Tm={c.pair.forward.tm.tm}°C")
        print(f"    REV {c.pair.reverse.sequence}  Tm={c.pair.reverse.tm.tm}°C")
        print(f"    ΔTm={c.pair.tm_difference}°C  heterodimer={c.heterodimer.dg:.2f} kcal/mol")
        for w in c.warnings:
            print(f"    ⚠ {w}")


def cmd_golden_gate(args):
    """Check Golden Gate overhangs."""
    result = check_golden_gate(
        args.overhangs,
        enzyme=args.enzyme,
        fidelity_set=args.fidelity_set,
    )

    print(f"\n{'=' * 60}")
    print(f"  GOLDEN GATE OVERHANG CHECK — {result.enzyme}")
    print(f"{'=' * 60}")

    for oh in result.overhangs:
        fid = '✓ NEB' if oh.is_in_fidelity_set else '  custom'
        pal = '⚠ PALINDROME' if oh.is_palindromic else ''
        print(f"  {oh.overhang}  rc={oh.reverse_complement}  GC={oh.gc_content:.0%}  {fid}  {pal}")

    print(f"\n  Unique: {'✓' if result.all_unique else '✗ COLLISION'}")

    for w in result.warnings:
        print(f"  ⚠ {w}")


def main():
    parser = argparse.ArgumentParser(
        prog='primerdesignr',
        description='Primer design and thermodynamic analysis',
    )
    parser.add_argument('--na', type=float, default=50.0, help='Na+ concentration (mM)')
    parser.add_argument('--mg', type=float, default=0.0, help='Mg2+ concentration (mM)')
    parser.add_argument('--dna', type=float, default=250.0, help='Primer concentration (nM)')

    sub = parser.add_subparsers(dest='command')

    # analyze
    p_analyze = sub.add_parser('analyze', help='Analyze primer sequences')
    p_analyze.add_argument('sequences', nargs='*', help='Primer sequences or NAME:SEQUENCE pairs')
    p_analyze.add_argument('--csv', help='Export results to CSV file')

    # pair
    p_pair = sub.add_parser('pair', help='Analyze a primer pair')
    p_pair.add_argument('forward', help='Forward primer sequence')
    p_pair.add_argument('reverse', help='Reverse primer sequence')

    # batch
    p_batch = sub.add_parser('batch', help='Analyze primers from a file')
    p_batch.add_argument('file', help='File with primer sequences')
    p_batch.add_argument('--csv', help='Export results to CSV file')

    # design
    p_design = sub.add_parser('design', help='Design PCR primer pairs from a template')
    p_design.add_argument('template', nargs='?', help='Template sequence; omit to read stdin')
    p_design.add_argument('--file', help='Read template sequence from a file')
    p_design.add_argument('--product-min', type=int, default=120, help='Minimum product size')
    p_design.add_argument('--product-max', type=int, default=500, help='Maximum product size')
    p_design.add_argument('--count', type=int, default=5, help='Number of candidate pairs')

    # golden-gate
    p_gg = sub.add_parser('golden-gate', help='Check Golden Gate overhangs')
    p_gg.add_argument('overhangs', nargs='+', help='4bp overhang sequences')
    p_gg.add_argument('--enzyme', default='BsaI', help='Type IIS enzyme (default: BsaI)')
    p_gg.add_argument('--fidelity-set', default='10_fragment', help='NEB fidelity set to check against')

    args = parser.parse_args()

    if args.command == 'analyze':
        cmd_analyze(args)
    elif args.command == 'pair':
        cmd_pair(args)
    elif args.command == 'batch':
        cmd_batch(args)
    elif args.command == 'design':
        cmd_design(args)
    elif args.command == 'golden-gate':
        cmd_golden_gate(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
