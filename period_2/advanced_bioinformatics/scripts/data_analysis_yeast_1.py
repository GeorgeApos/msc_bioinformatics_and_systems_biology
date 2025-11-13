#!/usr/bin/env python3
"""
Yeast IMW004 Variant Analysis Pipeline
"""

import subprocess
import sys
import os

REF = "files/yeast/chr3.fasta"
GFF = "files/yeast/chr3.gff"

# Paired-end reads
READ1 = "files/yeast/imw004-chr3_1.fastq"
READ2 = "files/yeast/imw004-chr3_2.fastq"

SAM = "files/yeast/imw004.sam"
BAM = "files/yeast/imw004.bam"
SORTED = "files/yeast/imw004.sorted.bam"
VCF = "files/yeast/variants.vcf"


def run(cmd):
    print(f"\n[RUN] {cmd}")
    ret = subprocess.run(cmd, shell=True)
    if ret.returncode != 0:
        sys.exit(f"Error running: {cmd}")


def mapping_and_variants():

    run(f"bwa index {REF}")

    # *** FIXED: use paired-end reads ***
    run(f"bwa mem {REF} {READ1} {READ2} > {SAM}")

    run(f"samtools view -S -b {SAM} > {BAM}")
    run(f"samtools sort {BAM} -o {SORTED}")
    run(f"samtools index {SORTED}")

    run(f"bcftools mpileup -Ou -f {REF} {SORTED} | bcftools call -mv -Ov -o {VCF}")

    print(f"\nVCF written to {VCF}\n")


def find_ady2_in_gff():
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            attributes = cols[8]
            if "ADY2" in attributes:
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6]
                print(f"ADY2 gene found at: {start}-{end} (strand {strand})")
                return start, end, strand

    sys.exit("ADY2 gene not found in GFF!")


def find_snps(vcf_file, gene_start, gene_end):
    print("\nChecking SNPs in ADY2...\n")
    snps = []

    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, _, _, _ = line.strip().split("\t", 7)
            pos = int(pos)

            if gene_start <= pos <= gene_end:
                snps.append((pos, ref, alt))
                print(f"SNP at {pos}: {ref} → {alt}")

    if not snps:
        print("No SNPs detected in ADY2.")
    else:
        print("\nTotal ADY2 SNPs found:", len(snps))

    return snps


def main():
    print("\n=== YEAST PIPELINE START ===")

    mapping_and_variants()

    gene_start, gene_end, strand = find_ady2_in_gff()
    find_snps(VCF, gene_start, gene_end)

    if strand == "-":
        print("\nNOTE: ADY2 is on the reverse strand.\n")

    print("\n=== PIPELINE FINISHED ===\n")


if __name__ == "__main__":
    main()
