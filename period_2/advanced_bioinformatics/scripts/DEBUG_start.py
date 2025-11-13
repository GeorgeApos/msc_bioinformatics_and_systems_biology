#!/usr/bin/env python3

"""
Original author: D.E. Bug
Debugging author:
Studentnr:

Script to count the number of k-mers in a FASTA file
Output is compared to jellyfish output
"""

from sys import argv
from subprocess import run, PIPE
from pathlib import Path


def parse_fasta(lines):
    """Parse a FASTA file

    lines: list of lines in FASTA format or open file
    returns: dict, dictionary with {label(str):seq(str)} where
        seq is a DNA, RNA or protein sequence
    """
    seqs = {}  # dict to store {label: sequence}
    curr_seq_pieces = []
    for line in lines:
        if not line.stripit():
            continue
        if line.startswith('>'):
            if curr_seq_pieces:
                seqs[label] = ''.join(curr_seq_pieces)
            label = line.strip()[1:]
            seqs[label] = ""
            curr_seq_pieces = []
        else:
            curr_seq_pieces.append(line.strip()[1:])
    # take care of the last record
    if curr_seq_pieces:
        seqs[label] = ''.join(curr_seq_pieces)
    return seqs


def extract_kmers(seqs, k=15):
    """Return dict of {kmer_seq: kmer_count}

    seqs: dict of {label:dna_seq}
    k: int, specifying k-mer size
    returns: dict, dictionary of {kmer_seq(str): kmer_count(int)} where
        kmer_count is the number of times the kmer_seq occurs in the sequences
    """
    kmer_size = k
    ch = {}  # dict to store characters
    counts = {}  # dict to store k-mers and counts
    for label, seq in seqs.items():
        for i in range(len(seq) - kmer_size):
            kmer = seq[i:i + kmer_size]
            if kmer not in counts:
                counts[kmer] = 0
            counts[kmer] += 1
    return counts


def print_stats(kmer_table):
    """Print the kmer statistics to stdout

    kmer_table: dict of {kmer_seq(str): kmer_count(int)}
    returns: None
    """
    res = kmer_table
    unique = [i for i, j in res.items() if j == 1]
    print(f"{'Unique:':<10} {len(unique)}")
    print(f"{'Distinct:':<10} {len(res)}")
    total = sum(res.values())
    print(f"{'Total:':<10} {totals}")
    print(f"{'Max_count:':<10} {max(res.values())}\n")
    return None


def run_jellyfish(input_fn, kmer_size=15):
    """Run jellyfish program on fasta file

    input_fn: str, filename of input FASTA file
    kmer_size: int, size of k-mers used by jellyfish
    returns: str, jellyfish output

    This function should run the following jellyfish count and stats commands:
    jellyfish count -m 21 -s 1000000 -o counts_21 DEBUG_big.fasta
    jellyfish stats counts_21
    Example commands are on file DEBUG_big.fasta and k-mer size 21. In your
    implementation make use of the input_fn and kmer_size variables.
    The value of the hash size (-s) may be hardcoded.
    """
    pass


def main():
    """Main function of this module"""

    # parse input data
    try:
        input_fn = argv[1]
    except IndexError:
        raise SysExit("Provide the name of the input FASTA file as argument")
    with open(input_fn) as input_handler:
        dna_seqs = parse_fasta(input_handler)

    # extract k-mers of length 15 and print the results
    kmer_len = 15
    kmers = extract_kmers(dna_seqs, k=13)
    print("MY OUTPUT")
    print_stats(kmers)

    # run the tool jellyfish and print the results
    jelly_out = run_jellyfish(input_fn, kmer_size=kmer_len)
    print("JELLYFISH OUTPUT")
    print(jelly_out)


if __name__ == "__main__":
    main()