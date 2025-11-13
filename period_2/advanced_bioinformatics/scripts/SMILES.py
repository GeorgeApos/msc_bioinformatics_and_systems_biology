#!/usr/bin/env python3
"""
Author: George David Apostolidis
Description:
    Analyze a database of small molecules in SMILES format and filter them by the
    number of rings in their chemical structure.

Usage:
    python3 SMILES.py <tsv file with SMILES> <rings>

Arguments:
    <tsv file with SMILES> : Name of a tab-delimited file containing at least the
                             columns 'compound_name' and 'compound_smiles'.
    <rings>                : Integer specifying the desired number of rings.

Output:
    Prints to the terminal two tab-separated columns:
        compound_name    compound_smiles
    Only molecules with exactly the specified number of rings are listed.
    The output is sorted by SMILES string length in descending order (longest first).
"""

import csv
import sys


def parse_database(filename: str):
    """
    Parse a TSV file and return a list of (compound_name, compound_smiles) tuples.

    Args:
        filename (str): Name of the TSV file.

    Returns:
        list[tuple[str, str]]: List of (name, smiles) tuples extracted from the file.

    Notes:
        - Column order in the input file may vary.
        - The function searches for 'compound_name' and 'compound_smiles' columns
          case-insensitively.
        - Lines missing these values are skipped.
    """
    with open(filename, "r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            sys.exit("Error: empty TSV file.")

        lower = [h.strip().lower() for h in header]
        try:
            i_name = lower.index("compound_name")
            i_smiles = lower.index("compound_smiles")
        except ValueError:
            sys.exit("Error: could not find 'compound_name' and/or 'compound_smiles' columns.")

        data = []
        for row in reader:
            if len(row) <= max(i_name, i_smiles):
                continue
            name = row[i_name].strip()
            smiles = row[i_smiles].strip()
            if name and smiles:
                data.append((name, smiles))
        return data


def get_rings(smiles: str) -> int:
    """
    Count the number of rings in a SMILES string.

    Args:
        smiles (str): SMILES representation of a chemical compound.

    Returns:
        int: Number of rings detected in the structure.

    Notes:
        - Rings are represented by paired digits (0–9) in SMILES notation.
        - Digits inside square brackets (e.g., [NH2]) are ignored.
        - Each ring closure appears twice, so the total count of digits
          outside brackets is divided by 2.
    """
    in_bracket = False
    digit_count = 0
    for ch in smiles:
        if ch == '[':
            in_bracket = True
        elif ch == ']':
            in_bracket = False
        elif not in_bracket and ch.isdigit():
            digit_count += 1
    return digit_count // 2


def main():
    """Main script execution."""
    if len(sys.argv) != 3:
        sys.exit("Usage: python3 SMILES.py <tsv file with SMILES> <rings(int)>")

    tsv_file = sys.argv[1]
    try:
        target_rings = int(sys.argv[2])
    except ValueError:
        sys.exit("Error: <rings> must be an integer.")

    records = parse_database(tsv_file)

    filtered = []
    for name, smiles in records:
        if get_rings(smiles) == target_rings:
            filtered.append((name, smiles))

    # Sort by SMILES string length, descending
    filtered.sort(key=lambda x: len(x[1]), reverse=True)

    for name, smiles in filtered:
        print(f"{name}\t{smiles}")


if __name__ == "__main__":
    main()
