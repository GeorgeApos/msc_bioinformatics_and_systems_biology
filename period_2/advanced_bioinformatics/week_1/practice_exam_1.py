#!/usr/bin/env python3
"""
Author: George David Apostolidiss
Description: Functions to calculate GC content and reverse complement of a DNA sequence
Usage: python3 my_exam.py velvet_15.fa chr3.fa
    dnasequence: str, DNA sequence over the alphabet A, C, G, T
"""

import sys
import subprocess

from jupyter_core.version import parts


def accept_args():
    """Accept command-line arguments.

    Returns:
        filenames: list of str, input filenames
    """

    if len(sys.argv) != 4:
        print("Usage: e.g. python3 my_exam.py velvet_15.fa chr3.fa")
        sys.exit(1)

    filenames = sys.argv[1:3]
    return filenames


def compare_assemblies(param, param1):

    """Compare two assembly files using diff.

    Args:
        param: str, first filename
        param1: str, second filename
    """

    cmd = [
        "lastz",
        "velvet_15.fa[multiple]",
        "chr3.fa[multiple]",
        "--chain",
        "--gapped",
        "--format=general:name1,start1,end1,name2,start2,end2,identity"
    ]

    with open("output.txt", "w") as outfile:
        subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)


def parse_output_to_extract_alignments(param, param1):

    alignments = []

    with open("output.txt", "r") as infile:
        for line in infile:
            parts = line.rstrip("\n").split("\t")
            if line.startswith("#") or parts[0].lower() == "name1":
                continue
            name1 = parts[0]
            start1 = int(parts[1])
            end1 = int(parts[2])
            name2 = parts[3]
            start2 = int(parts[4])
            end2 = int(parts[5])
            identity = parts[6]
            idPct = parts[7]

            alignments.append((name1, start1, end1, name2, start2, end2, identity, idPct))

    return alignments


def find_uncovered_regions(reference_genome, current_sequence):
    i, j = 0, 0
    uncovered_regions = {}
    count_bases = 0

    print("Uncovered regions:")

    while i < len(reference_genome) and j < len(current_sequence):
        if reference_genome[i] == current_sequence[j]:
            i += 1
            j += 1
        else:
            start = i
            char = reference_genome[i]
            count_bases += 1
            while i < len(reference_genome) and (j >= len(current_sequence) or reference_genome[i] != current_sequence[j]):
                i += 1
                char += reference_genome[i]
                count_bases += 1
            end = i - 1
            uncovered_regions[(start, end)] = char
            print(uncovered_regions[(start, end)] + " " + char)

    print("Number of uncovered regions:", len(uncovered_regions))
    print("Number of bases in uncovered regions:", count_bases)



def calculate_uncovered_regions(param, param1):

    reference_genome = ""

    with open(param1, "r") as infile:
        for line in infile:
            if not line.startswith(">"):
                reference_genome += line.rstrip("\n")

    with open(param, "r") as infile:
        current_sequence = ""

        for line in infile:
            if not line.startswith(">"):
                current_sequence += line.rstrip("\n")
            else:
                if current_sequence != "":
                    find_uncovered_regions(reference_genome, current_sequence)
                current_sequence = ""

def main():
    """Main function: entry point for the script."""
    filenames = accept_args()

    compare_assemblies(filenames[0], filenames[1])

    alignments = parse_output_to_extract_alignments(filenames[0], filenames[1])

    calculate_uncovered_regions(filenames[0], filenames[1])

if __name__ == "__main__":
    main()

