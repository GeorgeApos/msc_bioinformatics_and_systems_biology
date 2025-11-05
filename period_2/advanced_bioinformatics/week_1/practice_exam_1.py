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

    if len(sys.argv) != 3:
        print("Usage: e.g. python3 my_exam.py velvet_15.fa chr3.fa")
        sys.exit(1)

    filenames = sys.argv[1:]
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


def fasta_length(param):
    """Calculate the length of sequences in a FASTA file.

    Args:
        param: str, filename
    Returns:
        int: total length of sequences
    """

    total_length = 0
    with open(param, "r") as infile:
        for line in infile:
            if not line.startswith(">"):
                total_length += len(line.strip())
    return total_length


def read_lastz_intervals(path):
    """Read lastz output and extract intervals.

    Args:
        path: str, path to lastz output file
    Returns:
        list of tuples: (start, end) intervals
    """

    with open(path, "r") as infile:
        intervals = {}
        for line in infile:
            parts = line.rstrip("\n").split("\t")
            if line.startswith("#") or parts[0].lower() == "name1":
                continue
            ref = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            intervals.setdefault(ref, []).append((start, end))
        return intervals


def main():
    """Main function: entry point for the script."""
    filenames = accept_args()

    compare_assemblies(filenames[0], filenames[1])

    chr3_length = fasta_length(filenames[1])

    ivs = read_lastz_intervals("output.txt")
    print(len(ivs))

    pass


if __name__ == "__main__":
    main()