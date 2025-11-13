#!/usr/bin/env python3
"""
Author: George David Apostolidis

Description:
    This script parses a FASTQ file containing sequencing reads and identifies
    the sequence(s) with the greatest length. Each record in a FASTQ file
    consists of four lines:
        1. A header line beginning with '@' followed by the sequence label and metadata
        2. The nucleotide sequence (read)
        3. A '+' separator line
        4. A quality string, with one character per base indicating read quality

    The program reads all records, determines the maximum sequence length, and
    prints the label(s) of the read(s) that have this maximum length.

Usage:
    python3 FASTQ1.py <FILENAME.fastq>

Arguments:
    <FILENAME.fastq> : Path to the FASTQ file to be parsed.

Output:
    One or more sequence labels (without the '@' symbol) printed to standard output,
    corresponding to the longest sequence(s) in the file.

Example:
    $ python3 FASTQ1.py FASTQ1_sample.fastq
    1a427025-1eb6-40fe-9c8e-70ca779706e7
"""


from typing import List, Tuple
import sys


def parse_fastq(filename: str) -> List[Tuple[str, str]]:
    """
    Read a FASTQ file and return a list of (sequence, label) tuples.

    The "label" is the sequence identifier (the token immediately following '@'
    on the header line, up to the first whitespace). Any additional metadata
    on the header line is ignored.

    Parameters
    ----------
    filename : str
        Path to the FASTQ file.

    Returns
    -------
    List[Tuple[str, str]]
        A list of (sequence, label) tuples in file order.

    Raises
    ------
    ValueError
        If the file is malformed (e.g., a header not starting with '@' or a
        missing '+' line).
    """
    records: List[Tuple[str, str]] = []

    with open(filename, "r", encoding="utf-8") as fh:
        while True:
            header = fh.readline()
            if not header:  # EOF
                break

            if not header.startswith("@"):
                raise ValueError("Malformed FASTQ: header line does not start with '@'.")

            # The label is the first token after '@'
            label = header[1:].strip().split()[0]

            seq = fh.readline()
            if not seq:
                raise ValueError("Malformed FASTQ: missing sequence line.")
            seq = seq.strip()

            plus = fh.readline()
            if not plus:
                raise ValueError("Malformed FASTQ: missing '+' line.")
            if not plus.startswith("+"):
                raise ValueError("Malformed FASTQ: third line must start with '+'.")

            qual = fh.readline()
            if not qual:
                raise ValueError("Malformed FASTQ: missing quality line.")
            qual = qual.rstrip("\n")

            # Sanity check: length of quality string should equal sequence length.
            # We won't fail hard if it doesn't; just a gentle check.
            if len(qual) != len(seq):
                # You could raise here if the assignment requires strict validation:
                # raise ValueError("Quality length differs from sequence length.")
                pass

            records.append((seq, label))

    return records


def labels_longest_sequences(seq_label_list: List[Tuple[str, str]]) -> List[str]:
    """
    Given a list of (sequence, label) tuples, return all labels whose
    sequences have the maximum length.

    Parameters
    ----------
    seq_label_list : List[Tuple[str, str]]
        The parsed FASTQ records as (sequence, label).

    Returns
    -------
    List[str]
        Labels for all sequences that share the maximum length. Order is the
        same as input order.
    """
    if not seq_label_list:
        return []

    max_len = max(len(seq) for seq, _ in seq_label_list)
    return [label for seq, label in seq_label_list if len(seq) == max_len]


def main() -> None:
    """
    Entry point: parse arguments, read FASTQ, compute longest label(s),
    and print them (one per line).
    """
    if len(sys.argv) != 2:
        sys.exit(f"Usage: python3 {sys.argv[0]} <FILENAME.fastq>")

    filename = sys.argv[1]

    try:
        records = parse_fastq(filename)
    except (OSError, ValueError) as e:
        sys.exit(f"Error: {e}")

    longest = labels_longest_sequences(records)
    for label in longest:
        print(label)


if __name__ == "__main__":
    main()
