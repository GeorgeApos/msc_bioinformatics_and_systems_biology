#!/usr/bin/env python3
"""
Author: George David Apostolidis
Student ID: 1746227
Description: Comparing DNA sequences using Needleman-Wunsch alignment via EMBOSS Needle.
             The script takes two FASTA files containing DNA sequences and a gap penalty
                as input, performs global alignment using the Needleman-Wunsch algorithm,
                and outputs alignment metrics including identity percentage, Hamming distance,
                and alignment length.
Usage: Usage:python3 practice_exam_10_11_2025.py PTEref.fa PTErelated.fa <gap_penalty>
"""
import subprocess
import sys
from typing import List, Dict
from dataclasses import dataclass


@dataclass
class AlignmentMetrics:
    identity_percentage: float
    hamming_distance: int
    alignment_length: int

@dataclass
class AlignmentSummary:
    seq1_label: str
    seq2_label: str
    base_1: str
    base_2: str
    metrics: AlignmentMetrics

@dataclass
class NeedleOutput:
    summaries: List[AlignmentSummary]

def input_args():
    """
        Parse and validate command-line arguments.

        Parameters
        ----------
        None

        Returns
        -------
        filenames : List[str]
            List containing the two input FASTA filenames.
        gap_penalty : int
            Gap penalty for the Needleman-Wunsch alignment.
    """

    if len(sys.argv) != 4:
        print("Usage:python3 practice_exam_10_11_2025.py PTEref.fa PTErelated.fa <gap_penalty>")
        sys.exit(1)

    filenames = sys.argv[1:3]
    gap_penalty = int(sys.argv[3])
    return filenames, gap_penalty

def run_cmd(cmd: List[str]) -> subprocess.CompletedProcess:
    """
        Run a command using subprocess and handle errors.

        Parameters
        ----------
        cmd : List[str]
            Command and its arguments to be executed.

        Returns
        -------
        subprocess.CompletedProcess
            The result of the executed command.
    """

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd)}", file=sys.stderr)
        print(f"Return code: {e.returncode}", file=sys.stderr)
        print(f"Output: {e.output}", file=sys.stderr)
        sys.exit(1)


def parse_needle_output(OUTPUT_FILE, needle_output=None):
    """
        Parse the Needle output file to extract alignment summaries.

        Parameters
        ----------
        OUTPUT_FILE : str
            Path to the Needle output file.
        needle_output : NeedleOutput, optional
            Existing NeedleOutput object to populate (default is None).

        Returns
        -------
        NeedleOutput
            Parsed NeedleOutput object containing alignment summaries.
    """
    with open(OUTPUT_FILE, "r") as f:

        alignments = []
        label1=label2=None
        base_1=base_2=None

        for line in f:
            line = line.strip()
            if line.startswith("#======================================="):
                if label1 and label2 and base_1 and base_2:
                    alignments.append((label1, label2, base_1, base_2))
                label1=label2=None
                base_1=base_2=None
            elif line.startswith("# 1: "):
                label1 = line.split("# 1: ")[1].strip()
            elif line.startswith("# 2: "):
                label2 = line.split("# 2: ")[1].strip()
            elif not line.startswith("#") and line:
                parts = line.split()
                if len(parts) >= 3:
                    base_1 = parts[0]
                    base_2 = parts[2]

        if label1 and label2 and base_1 and base_2:
            alignments.append((label1, label2, base_1, base_2))
    summaries = []
    for label1, label2, base_1, base_2 in alignments:
        metrics = calculate_metrics(base_1, base_2)
        summary = AlignmentSummary(
            seq1_label=label1,
            seq2_label=label2,
            base_1=base_1,
            base_2=base_2,
            metrics=metrics
        )
        summaries.append(summary)
    needle_output = NeedleOutput(summaries=summaries)
    print()

    return needle_output


def calculate_metrics(aligned_seq1, aligned_seq2):
    """
        Calculate alignment metrics from two aligned sequences.

        Parameters
        ----------
        aligned_seq1 : str
            The first aligned sequence.
        aligned_seq2 : str
            The second aligned sequence.

        Returns
        -------
        AlignmentMetrics
            Calculated alignment metrics including identity percentage,
            Hamming distance, and alignment length.
    """
    alignment_length = len(aligned_seq1.replace('-', ''))
    identity_count = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
    hamming_distance = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a != b and a != '-' and b != '-')
    identity_percentage = (identity_count / alignment_length) * 100 if alignment_length > 0 else 0.0

    metrics = AlignmentMetrics(
        identity_percentage=identity_percentage,
        hamming_distance=hamming_distance,
        alignment_length=alignment_length
    )

    return metrics

def print_output(needle_output: NeedleOutput):
    """
        Print the alignment summaries from Needle output.

        Parameters
        ----------
        needle_output : NeedleOutput
            Parsed NeedleOutput object containing alignment summaries.

        Returns
        -------
        None
    """
    print("Sequence 1 Sequence 2 Length Hamm Ident")
    for summary in needle_output.summaries:
        print(f"{summary.seq1_label} {summary.seq2_label} "
              f"{summary.metrics.alignment_length} "
              f"{summary.metrics.hamming_distance} "
              f"{summary.metrics.identity_percentage:.2f}")

def main():
    """
        Main function to execute the Needleman-Wunsch alignment and process results.
    """
    OPEN_GAP_PENALTY = 8
    GAP_EXTENSION_PENALTY = 0.5
    OUTPUT_FILE = "out.needle"

    filenames, gap_penalty = input_args()

    if gap_penalty < 3 or gap_penalty > 10:
        gap_penalty = OPEN_GAP_PENALTY

    cmd = [
        "needle",
        "-asequence", filenames[0],
        "-bsequence", filenames[1],
        "-gapopen", str(gap_penalty),
        "-gapextend", str(GAP_EXTENSION_PENALTY),
        "-outfile", OUTPUT_FILE,
    ]

    run_cmd(cmd)

    needle_output = parse_needle_output(OUTPUT_FILE)

    print_output(needle_output)
    

if __name__ == "__main__":
    main()

