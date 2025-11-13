#!/usr/bin/env python3
"""
Author: George David Apostolidis (student no: XXXXXXX)
Description: Compare a Velvet assembly to a reference genome with LASTZ,
parse the alignments, and report UNCOVERED regions in the reference.

Usage:
    python3 my_exam.py velvet_15.fa chr3.fa
Output:
    Prints uncovered regions as "start:end sequence" (0-based, end inclusive),
    then the total number of regions and total number of uncovered bases.
"""

import os
import sys
import subprocess
from typing import List, Tuple, Iterable


LASTZ_OUT = "outlastz.txt"


def die(msg: str, code: int = 1) -> None:
    print(msg, file=sys.stderr)
    sys.exit(code)


def read_fasta_sequence(path: str) -> str:
    """Read a single- or multi-line FASTA into a single string (no headers)."""
    seq_chunks: List[str] = []
    try:
        with open(path, "r") as fh:
            for line in fh:
                if not line.startswith(">"):
                    seq_chunks.append(line.strip())
    except FileNotFoundError:
        die(f"Error: file not found: {path}")
    return "".join(seq_chunks)


def run_lastz(query_fa: str, ref_fa: str) -> None:
    """
    Run LASTZ once to produce outlastz.txt in 'general' format.
    If the output already exists, do not re-run (per exam guidance).
    """
    if os.path.exists(LASTZ_OUT):
        return

    # Ask for an explicit subset of 'general' fields so parsing is deterministic.
    # NOTE: LASTZ general coordinates are 1-based inclusive. We convert to 0-based half-open.
    fields = "name1,start1,end1,name2,start2,end2,identity,idPct"
    cmd = [
        "lastz",
        f"{query_fa}[multiple]",
        f"{ref_fa}[multiple]",
        "--chain",
        "--gapped",
        f"--format=general:{fields}",
    ]
    try:
        with open(LASTZ_OUT, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=True)
    except FileNotFoundError:
        die("Error: lastz not found in PATH on this system.")
    except subprocess.CalledProcessError as e:
        die(f"Error running lastz:\n{e.stderr or str(e)}")


def parse_lastz_general(path: str) -> Iterable[Tuple[int, int]]:
    """
    Parse LASTZ 'general' output and yield reference coverage intervals [start, end)
    in 0-based half-open coordinates, derived from fields start2,end2.

    We skip comments ('#...') and header line starting with 'name1'.
    """
    with open(path, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0].lower() == "name1":
                continue
            # columns: name1, start1, end1, name2, start2, end2, identity, idPct
            try:
                start2_1based = int(cols[4])
                end2_1based_inclusive = int(cols[5])
            except (IndexError, ValueError):
                continue  # skip malformed lines

            # Convert 1-based inclusive -> 0-based half-open:
            start = start2_1based - 1
            end = end2_1based_inclusive  # inclusive(1-based) -> exclusive(0-based)
            if start < end:
                yield (start, end)


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping/adjacent [start,end) intervals; return sorted, disjoint list."""
    sorted_iv = sorted(intervals)
    if not sorted_iv:
        return []
    merged: List[Tuple[int, int]] = [sorted_iv[0]]
    for s, e in sorted_iv[1:]:
        ps, pe = merged[-1]
        if s <= pe:  # overlap/adjacent
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def invert_intervals(covered: List[Tuple[int, int]], length: int) -> List[Tuple[int, int]]:
    """
    Given merged covered intervals and reference length, return uncovered intervals
    as 0-based half-open [start,end).
    """
    gaps: List[Tuple[int, int]] = []
    cursor = 0
    for s, e in covered:
        if cursor < s:
            gaps.append((cursor, s))
        cursor = max(cursor, e)
    if cursor < length:
        gaps.append((cursor, length))
    return gaps


def main() -> None:
    if len(sys.argv) != 3:
        die("Usage: python3 my_exam.py <assembly.fa> <reference.fa>")

    assembly_fa, reference_fa = sys.argv[1], sys.argv[2]

    # 1) Run LASTZ once
    run_lastz(assembly_fa, reference_fa)

    # 2) Read reference sequence (to slice uncovered regions later)
    ref_seq = read_fasta_sequence(reference_fa)
    ref_len = len(ref_seq)

    # 3) Parse coverage on reference and merge
    covered = merge_intervals(parse_lastz_general(LASTZ_OUT))

    # 4) Invert to get uncovered regions in 0-based half-open
    uncovered = invert_intervals(covered, ref_len)

    # 5) Report in exam format: "start:end sequence"
    #    The exam shows end as inclusive. Our intervals are half-open [s,e),
    #    so we print end-1 to match that style.
    print("Uncovered regions:")
    total_bases = 0
    for s, e in uncovered:
        seq = ref_seq[s:e]
        total_bases += (e - s)
        print(f"{s}:{e-1} {seq}")

    print()
    print(f"Number of uncovered regions: {len(uncovered)}")
    print(f"Number of uncovered bases: {total_bases}")


if __name__ == "__main__":
    main()
