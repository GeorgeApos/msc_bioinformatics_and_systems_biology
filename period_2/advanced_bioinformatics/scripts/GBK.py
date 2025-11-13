#!/usr/bin/env python3
"""
GBK.py

Parse a GenBank file (no external libs) containing one or more records and
produce:
  1) output.fasta  — sequences in FASTA, sorted by GC% (high → low)
     FASTA labels: "<accession> <organism>"
  2) output.tsv    — tab-delimited table with columns:
        Accession    Organism    GC(%) with two decimals (no % sign)    Length

Usage:
    python3 GBK.py <FILENAME.gb>
"""

import sys
from typing import List, Tuple


def gc_content(seq: str) -> float:
    """Compute GC percentage (0–100) for a nucleotide sequence.

    Counts only A/C/G/T (case-insensitive). Other symbols (e.g., N) are ignored
    in the denominator.

    Args:
        seq: Nucleotide sequence string.

    Returns:
        GC content as a float percentage in [0, 100].
    """
    if not seq:
        return 0.0
    s = seq.upper()
    gc = s.count('G') + s.count('C')
    atgc = sum(s.count(x) for x in ('A', 'C', 'G', 'T'))
    if atgc == 0:
        return 0.0
    return (gc / atgc) * 100.0


def wrap_fasta(seq: str, width: int = 60) -> List[str]:
    """Wrap a sequence string into fixed-width lines.

    Args:
        seq: Sequence to wrap.
        width: Maximum line width for the FASTA sequence.

    Returns:
        List of wrapped sequence lines (each up to `width` chars).
    """
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def parse_genbank(path: str) -> List[Tuple[str, str, str]]:
    """Parse a GenBank flatfile without external libraries.

    Extracts accession (first token on ACCESSION line), organism (from
    the ORGANISM line only, excluding lineage), and the sequence from the
    ORIGIN block of each record.

    Args:
        path: Path to the GenBank file (.gb/.gbk/.gbff) to parse.

    Returns:
        A list of (accession, organism, sequence) tuples for all records found.

    Raises:
        SystemExit: If the file is not found or permission is denied.
    """
    records: List[Tuple[str, str, str]] = []

    accession = None
    organism = None
    seq_chunks: List[str] = []
    in_sequence = False
    in_record = False

    try:
        with open(path, "r", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.rstrip("\n")

                # New record starts at LOCUS
                if line.startswith("LOCUS"):
                    if in_record and accession and organism:
                        seq = "".join(seq_chunks).replace(" ", "").upper()
                        records.append((accession, organism, seq))
                    accession = None
                    organism = None
                    seq_chunks = []
                    in_sequence = False
                    in_record = True
                    continue

                if not in_record:
                    continue

                if line.startswith("ACCESSION"):
                    parts = line.split()
                    if len(parts) >= 2:
                        accession = parts[1]
                    continue

                if line.startswith("  ORGANISM"):
                    try:
                        organism = line.split("ORGANISM", 1)[1].strip()
                    except Exception:
                        organism = ""
                    continue

                if line.startswith("ORIGIN"):
                    in_sequence = True
                    continue

                if line.strip() == "//":
                    seq = "".join(seq_chunks).replace(" ", "").upper()
                    seq = "".join(ch for ch in seq if ch.isalpha())
                    if accession and organism:
                        records.append((accession, organism, seq))
                    accession = None
                    organism = None
                    seq_chunks = []
                    in_sequence = False
                    in_record = False
                    continue

                if in_sequence:
                    cleaned = "".join(ch for ch in line if ch.isalpha())
                    if cleaned:
                        seq_chunks.append(cleaned)

        if in_record and accession and organism:
            seq = "".join(seq_chunks).replace(" ", "").upper()
            seq = "".join(ch for ch in seq if ch.isalpha())
            records.append((accession, organism, seq))

    except FileNotFoundError:
        sys.exit(f"Error: file not found: {path}")
    except PermissionError:
        sys.exit(f"Error: permission denied: {path}")

    return records


def write_outputs(records: List[Tuple[str, str, str]]) -> None:
    """Compute stats, sort by GC%, and write output files.

    Produces:
      - output.fasta: FASTA entries sorted by GC% (high → low),
        header as ">accession organism", wrapped at 60 chars.
      - output.tsv:   Tab-delimited rows: accession, organism, GC (two decimals, no %), length.

    Args:
        records: List of (accession, organism, sequence) tuples.

    Returns:
        None. Writes `output.fasta` and `output.tsv` in the current directory.
    """
    stats = []
    for acc, org, seq in records:
        length = len(seq)
        gc = gc_content(seq)
        stats.append((acc, org, seq, gc, length))

    stats.sort(key=lambda x: x[3], reverse=True)

    with open("output.fasta", "w", encoding="utf-8") as ffa:
        for acc, org, seq, gc, length in stats:
            ffa.write(f">{acc} {org}\n")
            for line in wrap_fasta(seq, 60):
                ffa.write(line + "\n")

    with open("output.tsv", "w", encoding="utf-8") as ftsv:
        for acc, org, seq, gc, length in stats:
            ftsv.write(f"{acc}\t{org}\t{gc:.2f}\t{length}\n")


def main() -> None:
    """Entry point for the GBK → FASTA/TSV converter.

    Reads the GenBank file path from argv, parses records, and writes outputs.
    """
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 GBK.py <FILENAME.gb>")
    infile = sys.argv[1]
    records = parse_genbank(infile)
    if not records:
        sys.exit("Error: no valid records parsed from GenBank file.")
    write_outputs(records)


if __name__ == "__main__":
    main()
