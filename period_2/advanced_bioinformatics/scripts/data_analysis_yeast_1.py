#!/usr/bin/env python3

"""
Author: George Apostolidis
Description: Pipeline for mapping IMW004 reads, calling variants and locating ADY2 SNP.
Usage: python3 pipeline.py
"""

from sys import argv
import subprocess
import os.path

# File paths (adjust if needed)
REF = "files/yeast/chr3.fasta"
GFF = "files/yeast/chr3.gff"

READ1 = "files/yeast/imw004-chr3_1.fastq"
READ2 = "files/yeast/imw004-chr3_2.fastq"

SAM = "files/yeast/imw004.sam"
BAM = "files/yeast/imw004.bam"
SORTED = "files/yeast/imw004.sorted.bam"
VCF = "files/yeast/variants.vcf"

# Path to VarScan (depends on your local installation)
VARSCAN = "/mnt/local_scratch/BIF30806/prog/varscan/VarScan.v2.3.9.jar"


def find_ady2_in_gff(gff_file):
    """Return the genomic region of ADY2 (start, end, strand)."""
    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            if "ADY2" in cols[8]:
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6]
                return start, end, strand
    return None


def find_snps(vcf_file, start, end):
    """Return all SNPs inside the ADY2 gene region."""
    hits = []
    with open(vcf_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            if start <= pos <= end:
                hits.append((pos, ref, alt))
    return hits


def main():

    commands = []

    # ---- 1. Index reference genome ----
    commands.append(f"bwa index {REF}")

    # ---- 2. Map IMW004 reads ----
    commands.append(f"bwa mem {REF} {READ1} {READ2} > {SAM}")

    # ---- 3. Convert SAM → BAM ----
    commands.append(f"samtools view -S -b {SAM} > {BAM}")

    # ---- 4. Sort BAM ----
    commands.append(f"samtools sort {BAM} -o {SORTED}")

    # ---- 5. Index BAM ----
    commands.append(f"samtools index {SORTED}")

    # ---- 6. Call variants using bcftools ----
    commands.append(
        f"bcftools mpileup -Ou -f {REF} {SORTED} | "
        f"bcftools call -mv -Ov -o {VCF}"
    )

    # ---- Show planned commands ----
    print("Getting ready to execute the following commands:\n")
    for c in commands:
        print(c)

    # ---- Execute commands ----
    for cmd in commands:
        print("\nExecuting:", cmd)
        subprocess.run(cmd, shell=True, check=True)
        print("Success")

    # ---- After running the pipeline: parse GFF and VCF ----
    print("\n=== ANALYSIS OF ADY2 ===")

    region = find_ady2_in_gff(GFF)
    if region is None:
        print("ADY2 not found in GFF!")
        return
    start, end, strand = region
    print(f"ADY2 gene region: {start}-{end} (strand {strand})")

    snps = find_snps(VCF, start, end)
    if not snps:
        print("No SNPs found in ADY2.")
    else:
        print("SNPs detected in ADY2:")
        for pos, ref, alt in snps:
            print(f" - Position {pos}: {ref} -> {alt}")

    print("\nPipeline finished.\n")


if __name__ == "__main__":
    main()
