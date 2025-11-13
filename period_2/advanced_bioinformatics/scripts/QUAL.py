#!/usr/bin/env python3
"""
QUAL: Per-base quality assessment and trimming for FASTQ files.

Performs:
  1) FASTQ parsing and quality score conversion (Phred+33 or Phred+64).
  2) Read length statistics and per-base average quality.
  3) External trimming with fastq_quality_trimmer (-t 30, auto -Q flag).
  4) Post-trim statistics and comparison table.

Usage:
    python3 QUAL.py <input.fq>
"""

from __future__ import annotations
from pathlib import Path
from typing import Iterable, List, Tuple
import shutil
import subprocess
import sys
import os


# ------------------------ FASTQ parsing & quality ------------------------ #

def parse_fastq(filepath: Path) -> List[Tuple[str, str]]:
    """Parse a FASTQ file and return (sequence, quality) pairs.

    Args:
        filepath (Path): Path to the FASTQ file to be parsed.

    Returns:
        List[Tuple[str, str]]: A list of (sequence, quality) tuples.

    Raises:
        ValueError: If FASTQ structure is malformed.
    """
    out: List[Tuple[str, str]] = []
    with filepath.open("r", encoding="utf-8") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            if not header.startswith("@"):
                raise ValueError("Malformed FASTQ: header line must start with '@'.")

            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()

            if not (seq and plus and qual):
                raise ValueError("Malformed FASTQ: incomplete record.")
            if not plus.startswith("+"):
                raise ValueError("Malformed FASTQ: third line must start with '+'.")

            out.append((seq.strip(), qual.strip()))
    return out


def detect_quality_offset(qualities: List[str]) -> int:
    """Detect Phred offset (33 or 64) from quality lines.

    Args:
        qualities (List[str]): List of raw quality strings from a FASTQ file.

    Returns:
        int: Detected Phred offset (33 or 64).
    """
    if not qualities:
        return 33
    min_code = min(ord(ch) for q in qualities for ch in q)
    return 64 if min_code >= 64 else 33


def convert_quality_line(quality_string: str, q_offset: int = 64) -> List[int]:
    """Convert a quality string to integer Phred scores.

    Args:
        quality_string (str): ASCII quality string from FASTQ.
        q_offset (int, optional): ASCII offset (33 or 64). Defaults to 64.

    Returns:
        List[int]: List of integer quality scores, clamped to [0, 41].
    """
    scores = []
    for ch in quality_string:
        q = ord(ch) - q_offset
        if q < 0:
            q = 0
        elif q > 41:
            q = 41
        scores.append(q)
    return scores


# ----------------------------- Stats helpers ----------------------------- #

def length_stats(seqs: Iterable[str]) -> Tuple[int, int, float]:
    """Compute minimum, maximum, and average read length.

    Args:
        seqs (Iterable[str]): Collection of sequence strings.

    Returns:
        Tuple[int, int, float]: (min_length, max_length, average_length).
    """
    lengths = [len(s) for s in seqs]
    if not lengths:
        raise ValueError("No sequences found.")
    return min(lengths), max(lengths), sum(lengths) / len(lengths)


def calculate_avg_q_per_position(qualities: List[List[int]]) -> List[float]:
    """Compute average quality at each base position.

    Args:
        qualities (List[List[int]]): List of lists of per-base quality scores.

    Returns:
        List[float]: Average quality per position.
    """
    if not qualities:
        return []
    max_len = max(len(q) for q in qualities)
    sums = [0] * max_len
    counts = [0] * max_len
    for q in qualities:
        for i, v in enumerate(q):
            sums[i] += v
            counts[i] += 1
    return [sums[i] / counts[i] if counts[i] else 0.0 for i in range(max_len)]


# -------------------------- Trimming integration -------------------------- #

def run_fastq_quality_trimmer(input_fq: Path, output_fq: Path, threshold: int, q_offset: int) -> None:
    """Run fastq_quality_trimmer to trim low-quality tails.

    Args:
        input_fq (Path): Input FASTQ file.
        output_fq (Path): Output trimmed FASTQ file (e.g., 'trimmed.fq').
        threshold (int): Quality threshold for trimming.
        q_offset (int): Phred offset (33 or 64) for encoding.

    Raises:
        SystemExit: If the tool fails or output file not created.
    """
    if shutil.which("fastq_quality_trimmer") is None:
        sys.exit("Error: fastq_quality_trimmer not found in PATH")

    cmd = [
        "fastq_quality_trimmer",
        "-t", str(threshold),
        "-Q", str(q_offset),
        "-i", str(input_fq),
        "-o", str(output_fq),
    ]
    try:
        res = subprocess.run(cmd, check=False, capture_output=True, text=True)
    except Exception as e:
        sys.exit(f"Error running fastq_quality_trimmer: {e}")

    if res.returncode != 0:
        msg = res.stderr.strip() or "Unknown error from fastq_quality_trimmer"
        sys.exit(f"Error running fastq_quality_trimmer: {msg}")

    if (not output_fq.is_file()) or os.path.getsize(output_fq) == 0:
        sys.exit("Error: trimmed.fq not created or empty")


# ---------------------------------- Main ---------------------------------- #

def main() -> None:
    """Main function: computes stats before and after trimming."""
    if len(sys.argv) != 2:
        sys.exit(f"Usage: python3 {Path(sys.argv[0]).name} <INPUT.fq>")

    input_path = Path(sys.argv[1]).resolve()
    if not input_path.is_file():
        sys.exit("Error: input file not found")

    records = parse_fastq(input_path)
    seqs = [s for s, _ in records]
    quals_raw = [q for _, q in records]

    q_offset = detect_quality_offset(quals_raw)
    quals_converted = [convert_quality_line(q, q_offset) for q in quals_raw]

    orig_min, orig_max, orig_avg = length_stats(seqs)
    orig_avgs = calculate_avg_q_per_position(quals_converted)

    trimmed_path = Path.cwd() / "trimmed.fq"
    if not trimmed_path.exists() or os.path.getsize(trimmed_path) == 0:
        run_fastq_quality_trimmer(input_path, trimmed_path, threshold=30, q_offset=q_offset)

    trimmed_records = parse_fastq(trimmed_path)
    trimmed_seqs = [s for s, _ in trimmed_records]
    trimmed_quals_raw = [q for _, q in trimmed_records]
    trimmed_quals_converted = [convert_quality_line(q, q_offset) for q in trimmed_quals_raw]

    trim_min, trim_max, trim_avg = length_stats(trimmed_seqs)
    trim_avgs = calculate_avg_q_per_position(trimmed_quals_converted)

    print(f"ORIGINAL: min={orig_min}, max={orig_max}, avg={orig_avg:.2f}")
    print(f"TRIMMED: min={trim_min}, max={trim_max}, avg={trim_avg:.2f}")

    max_positions = max(len(orig_avgs), len(trim_avgs))
    for pos in range(max_positions):
        o = orig_avgs[pos] if pos < len(orig_avgs) else 0.0
        t = trim_avgs[pos] if pos < len(trim_avgs) else 0.0
        imp = t - o
        print(f"{pos+1}\t{o:.2f}\t{t:.2f}\t{imp:.2f}")


if __name__ == "__main__":
    main()
