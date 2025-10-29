#!/usr/bin/env python3
"""
Author: George David Apostolidis
Description: Find and translate all protein sequences encoded by open reading frames (ORFs)
             in a DNA string. The script scans all six reading frames (3 forward and
             3 from the reverse complement) and prints each resulting protein sequence
             on a separate line. Uses the standard genetic code with ATG as the only
             start codon and TAA, TAG, TGA as stop codons.
Usage: python3 ORFs.py [DNA string]
"""

import sys

def genetic_code_rna():
    """Generate the standard genetic code (RNA codons).

    returns: dict, mapping of RNA codons (e.g., 'AUG') to amino acids (single-letter).
    """
    first  = 16 * 'U' + 16 * 'C' + 16 * 'A' + 16 * 'G'
    second = (4 * 'U' + 4 * 'C' + 4 * 'A' + 4 * 'G') * 4
    third  = 16 * 'UCAG'
    prot   = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    return {first[i] + second[i] + third[i]: prot[i] for i in range(len(first))}

CODON_TABLE_RNA = genetic_code_rna()
START_DNA = "ATG"
STOP_DNA = {"TAA", "TAG", "TGA"}

def to_rna(dna):
    """Convert a DNA string to RNA by replacing T with U.

    dna: str, DNA sequence over the alphabet A, C, G, T
    returns: str, RNA sequence over A, C, G, U
    """
    return dna.upper().replace("T", "U")

def reverse_complement(dna):
    """Create the reverse complement of a DNA sequence.

    dna: str, DNA sequence over the alphabet A, C, G, T
    returns: str, reverse complement DNA sequence (A<->T, C<->G), reversed order
    """
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    dna = dna.upper()
    return "".join(comp[b] for b in reversed(dna))

def translate_dna(dna_orf):
    """Translate a DNA ORF into a protein sequence (without terminal stop).

    dna_orf: str, DNA sequence starting at ATG and ending at a stop codon (inclusive)
    returns: str, amino acid sequence (single-letter), excluding the stop symbol
    """
    rna = to_rna(dna_orf)
    aa = []
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        aa_sym = CODON_TABLE_RNA.get(codon)
        if aa_sym is None:
            # Should not happen for valid ORFs; skip gracefully
            break
        if aa_sym == "*":  # stop codon at the end of ORF
            break
        aa.append(aa_sym)
    return "".join(aa)

def find_orfs_in_frame(dna, frame):
    """Find all ORFs in a single forward reading frame and return their protein translations.

    dna: str, DNA sequence (A, C, G, T)
    frame: int, reading frame offset (0, 1, or 2)
    returns: list[str], protein sequences translated from each ORF found
    """
    proteins = []
    n = len(dna)
    i = frame
    while i <= n - 3:
        codon = dna[i:i+3]
        if codon == START_DNA:
            # Found a start; search for the next in-frame stop
            j = i + 3
            while j <= n - 3:
                next_codon = dna[j:j+3]
                if next_codon in STOP_DNA:
                    orf_dna = dna[i:j+3]  # inclusive stop
                    proteins.append(translate_dna(orf_dna))
                    break
                j += 3
            # Continue scanning allowing overlapping ORFs
        i += 3
    return proteins

def find_all_orfs(dna):
    """Find all ORFs in all six reading frames and return their protein translations.

    dna: str, DNA sequence (A, C, G, T)
    returns: list[str], protein sequences from all ORFs (duplicates allowed)
    """
    dna = dna.upper()
    rev = reverse_complement(dna)
    all_proteins = []

    # Forward frames: 0, 1, 2
    for frame in (0, 1, 2):
        all_proteins.extend(find_orfs_in_frame(dna, frame))

    # Reverse-complement frames: 0, 1, 2
    for frame in (0, 1, 2):
        all_proteins.extend(find_orfs_in_frame(rev, frame))

    return all_proteins

def main():
    """Main function to parse input, find ORFs, and print translations.

    Reads a single DNA string from the command line and prints each translated
    protein (one per line). Prints nothing else.
    """
    if len(sys.argv) != 2:
        # Per spec, exactly one DNA string argument
        sys.exit(0)

    dna = sys.argv[1].strip().upper()
    # (Optional robustness) keep only valid DNA symbols
    # dna = "".join(b for b in dna if b in "ACGT")

    proteins = find_all_orfs(dna)
    for p in proteins:
        if p:  # print even if duplicates; skip empty strings just in case
            print(p)

if __name__ == "__main__":
    main()
