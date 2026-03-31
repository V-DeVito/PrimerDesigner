#!/usr/bin/env python3
"""
Fetch ssNA sequences from RCSB PDB for Binet et al. 2023 benchmark.

Run on your local machine (needs internet access):
    cd primerdesignr
    python3 tests/fetch_binet_sequences.py

Reads:   tests/data/binet2023_full.csv        (538 entries, no sequences)
Writes:  tests/data/binet2023_full_seqs.csv   (same + 'sequence' column)

Requires: requests  (pip install requests)
"""

import csv
import time
import sys

try:
    import requests
except ImportError:
    print("Install requests:  pip install requests")
    sys.exit(1)


def fetch_sequence_from_pdb(pdb_id, na_type, expected_size):
    """
    Fetch the nucleic acid sequence from a PDB entry.

    Tries polymer entities 1-4, picks the one whose type matches
    (polydeoxyribonucleotide for DNA, polyribonucleotide for RNA)
    and whose length matches expected_size.
    """
    na_keyword = "polydeoxyribonucleotide" if na_type == "DNA" else "polyribonucleotide"

    best_match = (None, None)

    for entity_id in range(1, 5):
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
        try:
            r = requests.get(url, timeout=10)
            if r.status_code != 200:
                continue
            data = r.json()
            entity_type = data.get("entity_poly", {}).get("type", "")
            seq = data.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
            seq = seq.replace("\n", "")

            if not seq:
                continue

            # Exact match: right NA type AND right length
            if na_keyword in entity_type.lower() and len(seq) == expected_size:
                return seq, entity_type

            # Partial match: right NA type, wrong length (save as fallback)
            if na_keyword in entity_type.lower() and best_match[0] is None:
                best_match = (seq, entity_type)

            # Length match but type unclear
            if len(seq) == expected_size and best_match[0] is None:
                best_match = (seq, entity_type)

        except Exception:
            continue

    return best_match


def main():
    input_file = "tests/data/binet2023_full.csv"
    output_file = "tests/data/binet2023_full_seqs.csv"

    with open(input_file) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    print(f"Fetching sequences for {len(rows)} PDB entries...")
    print(f"  DNA: {sum(1 for r in rows if r['na_type'] == 'DNA')}")
    print(f"  RNA: {sum(1 for r in rows if r['na_type'] == 'RNA')}")
    print()

    success = 0
    check = 0
    failed = 0

    for i, row in enumerate(rows):
        pdb_id = row["pdb"]
        na_type = row["na_type"]
        expected_size = int(row["size"])

        seq, entity_type = fetch_sequence_from_pdb(pdb_id, na_type, expected_size)

        if seq and len(seq) == expected_size:
            row["sequence"] = seq.upper()
            print(f"  ✓ {pdb_id}: {len(seq)} nt {na_type} ({entity_type})")
            success += 1
        elif seq:
            row["sequence"] = f"CHECK:{seq[:40]}..."
            print(f"  ? {pdb_id}: got {len(seq)} nt, expected {expected_size} — check manually")
            check += 1
        else:
            row["sequence"] = ""
            print(f"  ✗ {pdb_id}: fetch failed")
            failed += 1

        # Progress every 50
        if (i + 1) % 50 == 0:
            print(f"  --- {i + 1}/{len(rows)} done ---")

        time.sleep(0.25)

    # Write output
    fieldnames = list(rows[0].keys())
    if "sequence" not in fieldnames:
        fieldnames.append("sequence")

    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n{'='*50}")
    print(f"Done: {success} fetched, {check} need review, {failed} failed")
    print(f"Output: {output_file}")

    if check > 0:
        print(f"\nTo fix CHECK entries:")
        print(f"  1. Go to rcsb.org/structure/XXXX")
        print(f"  2. Find the nucleic acid chain")
        print(f"  3. Copy the sequence into the CSV")

    print(f"\nNext: bring {output_file} back to Claude for the benchmark run")


if __name__ == "__main__":
    main()
