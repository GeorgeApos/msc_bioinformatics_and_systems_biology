#!/usr/bin/env python3
"""
Author: George David Apostolidiss
Description: Functions to calculate GC content and reverse complement of a DNA sequence
Usage: python3 DNA.py dnasequence
    dnasequence: str, DNA sequence over the alphabet A, C, G, T
"""

# import statements
import sys


def calculate_gc_content(dnaseq):
    """Calculate the GC percentage in a DNA sequence

    dnaseq: str, DNA sequence over the alphabet A, C, G, T
    returns: float, percentage of GC in a DNA sequence
    """

    count = 0
    for base in dnaseq:
        if base == 'G' or base == 'C':
            count += 1
    gc_content = (count / len(dnaseq)) * 100
    return gc_content

    # Alternative implementation using str.count()
    # gc_count = dnaseq.count('G') + dnaseq.count('C')
    # gc_content = (gc_count / len(dnaseq)) * 100
    # return gc_content

def reverse_complement(dnaseq):
    """Create the reverse complement of a DNA sequence

    dnaseq: str, DNA sequence over the alphabet A, C, G, T
    returns: str, reverse complement of the DNA sequence
    """

    # Implementation using a for loop
    # rev_comp = ''
    # for base in reversed(dnaseq):
    #     if base == 'A':
    #         rev_comp += 'T'
    #     elif base == 'T':
    #         rev_comp += 'A'
    #     elif base == 'C':
    #         rev_comp += 'G'
    #     elif base == 'G':
    #         rev_comp += 'C'
    # return rev_comp

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_comp = ''.join(complement[base] for base in reversed(dnaseq))
    return rev_comp

    # Alternative implementation using str.translate()
    # translation_table = str.maketrans('ACGT', 'TGCA')
    # rev_comp = dnaseq.translate(translation_table)[::-1]
    # return rev_comp

def main():
    """Main function"""

    # Command Line Argument Checking
    if len(sys.argv) != 2:
        print("Usage: python3 DNA.py dnasequence")
        sys.exit(1)

    # Validate DNA sequence
    for base in sys.argv[1]:
        if base not in "ACGT":
            print("Error: DNA sequence must contain only A, C, G, T")
            sys.exit(1)

    # Get DNA sequence from command line argument
    sequence = sys.argv[1]

    # Call functions and print results
    print(calculate_gc_content(sequence))
    print(reverse_complement(sequence))

if __name__ == "__main__":
    main()