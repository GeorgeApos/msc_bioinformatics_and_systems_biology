#!/usr/bin/env python3
"""
Original author: D.E. Bug
Debugging author: <your name>
Studentnr: <your student number>

Script to count the number of k-mers in a FASTA file
Output is compared to jellyfish output
"""

from pathlib import Path
from subprocess import run, PIPE
import sys


def parse_fasta(lines):
    """Parse a FASTA file.

    Args:
        lines: iterable of lines (an open file handle or list of strings)

    Returns:
        dict: {label(str): seq(str)}
    """
    seqs = {}
    label = None
    curr_seq_pieces = []

    for raw in lines:
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            # flush previous
            if label is not None:
                seqs[label] = ''.join(curr_seq_pieces)
            # new label (use whole header minus '>' or just first token—either is fine;
            # jellyfish doesn’t depend on labels)
            label = line[1:]
            curr_seq_pieces = []
        else:
            curr_seq_pieces.append(line)

    # flush last record
    if label is not None:
        seqs[label] = ''.join(curr_seq_pieces)

    return seqs


def extract_kmers(seqs, k=15):
    """Return dict of {kmer_seq: kmer_count}.

    Args:
        seqs: dict {label: dna_seq}
        k: k-mer size

    Returns:
        dict: {kmer(str): count(int)}
    """
    counts = {}
    for _, seq in seqs.items():
        n = len(seq)
        if n < k:
            continue
        # inclusive end: n - k + 1 windows
        for i in range(n - k + 1):
            kmer = seq[i:i + k]
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def print_stats(kmer_table):
    """Print the k-mer statistics to stdout.

    Args:
        kmer_table: dict {kmer_seq: kmer_count}
    """
    unique_count = sum(1 for c in kmer_table.values() if c == 1)
    distinct = len(kmer_table)
    total = sum(kmer_table.values())
    max_count = max(kmer_table.values()) if kmer_table else 0

    print(f"{'Unique:':<10} {unique_count}")
    print(f"{'Distinct:':<10} {distinct}")
    print(f"{'Total:':<10} {total}")
    print(f"{'Max_count:':<10} {max_count}\n")


def run_jellyfish(input_fn, kmer_size=15):
    """Run jellyfish (count + stats) and return stats output as a string.

    Commands executed:
        jellyfish count -m <k> -s 1000000 -o counts_<k> <input_fn>
        jellyfish stats counts_<k>

    Args:
        input_fn: path to FASTA
        kmer_size: k-mer size

    Returns:
        str: the stdout of `jellyfish stats`
    """
    out_fn = f"counts_{kmer_size}"

    # Count (suppress normal output; capture stderr too to avoid clutter)
    run(
        ["jellyfish", "count",
         "-m", str(kmer_size),
         "-s", "1000000",
         "-o", out_fn,
         str(input_fn)],
        check=True,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )

    # Stats (capture stdout to print later)
    stats_proc = run(
        ["jellyfish", "stats", out_fn],
        check=True,
        stdout=PIPE,
        stderr=PIPE,
        text=True,
    )
    return stats_proc.stdout.rstrip("\n")


def main():
    """Main entry point."""
    try:
        input_fn = sys.argv[1]
    except IndexError:
        sys.exit("Provide the name of the input FASTA file as argument")

    # parse input data
    with open(input_fn) as input_handler:
        dna_seqs = parse_fasta(input_handler)

    # extract k-mers of length 15 and print the results
    kmer_len = 15
    kmers = extract_kmers(dna_seqs, k=kmer_len)

    print("MY OUTPUT")
    print_stats(kmers)

    # run jellyfish and print the results
    jelly_out = run_jellyfish(input_fn, kmer_size=kmer_len)
    print("JELLYFISH OUTPUT")
    print(jelly_out)


if __name__ == "__main__":
    main()
