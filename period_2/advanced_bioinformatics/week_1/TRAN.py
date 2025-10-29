#!/usr/bin/env python3
"""
Author: George David Apostolidis
Description: Script for translating RNA sequences into protein sequences using the standard genetic code.
Usage: python3 TRAN.py [RNA string list]
    Each RNA string must:
      - Start with a start codon (AUG)
      - End with a stop codon (UAA, UAG, or UGA)
      - Contain no internal stop codons
      - Have a length divisible by 3
    The script prints one translated protein sequence per valid RNA string.
"""

import sys

def genetic_code():
    """Generate the standard genetic code as a dictionary.

    returns: dict, mapping of RNA codons (e.g., 'AUG') to amino acids (single-letter code).
    """
    first  = 16 * 'U' + 16 * 'C' + 16 * 'A' + 16 * 'G'
    second = (4 * 'U' + 4 * 'C' + 4 * 'A' + 4 * 'G') * 4
    third  = 16 * 'UCAG'
    prot   = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    return {first[i] + second[i] + third[i]: prot[i] for i in range(len(first))}


CODON_TABLE = genetic_code()
STOP_CODONS = {'UAA', 'UAG', 'UGA'}


def is_valid(rna: str) -> bool:
    """Check if an RNA sequence is a valid coding sequence (CDS).

    rna: str, RNA sequence (A, U, G, C)
    returns: bool, True if the sequence is valid, False otherwise.

    A valid CDS satisfies:
        - Length is divisible by 3
        - Starts with AUG
        - Ends with a stop codon (UAA, UAG, or UGA)
        - Contains only valid codons found in the standard genetic code
        - Has no internal stop codons
    """
    rna = rna.strip().upper()
    if len(rna) == 0 or (len(rna) % 3) != 0:
        return False
    if rna[:3] != 'AUG':
        return False
    if rna[-3:] not in STOP_CODONS:
        return False

    # Check all codons valid and no internal stop codons
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        if codon not in CODON_TABLE:
            return False
        if i not in (0, len(rna) - 3) and CODON_TABLE[codon] == '*':
            return False
    return True


def translate(rna: str) -> str:
    """Translate a valid RNA sequence into a protein sequence.

    rna: str, valid RNA sequence (A, U, G, C)
    returns: str, protein sequence (amino acid single-letter codes), excluding the terminal stop codon.
    """
    rna = rna.strip().upper()
    aa = []
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        aa_sym = CODON_TABLE[codon]
        if aa_sym == '*':  # terminal stop at the end
            break
        aa.append(aa_sym)
    return ''.join(aa)


def main():
    """Main program function.

    Reads RNA sequences from the command line, validates them,
    translates valid sequences, and prints each resulting protein on a separate line.
    """
    if len(sys.argv) < 2:
        print("Usage: python3 TRAN.py [RNA string list]")
        sys.exit(1)

    for rna in sys.argv[1:]:
        rna = rna.strip().upper()
        if is_valid(rna):
            print(translate(rna))


if __name__ == '__main__':
    main()
