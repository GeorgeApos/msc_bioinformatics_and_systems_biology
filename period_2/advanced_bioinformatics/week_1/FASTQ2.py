#!/usr/bin/env python3
"""
Author: George David Apostolidis

FASTQ2
Cross-platform version.

Description:
    - If input is .fastq: parse and print the label(s) of the longest sequence(s)
    - If input is .zip: unzip using Linux 'unzip' (quietly) or Python's zipfile on Windows,
      extract into the same folder as this script, and then parse the resulting .fastq

Usage:
    python3 FASTQ2.py <path to FILENAME.fastq>
    python3 FASTQ2.py <path to FILENAME.zip>

If the argument is not a valid file, the program exits with:
    Error: invalid file
"""

from __future__ import annotations
from pathlib import Path
from typing import List, Tuple
import subprocess
import platform
import zipfile
import sys


# ----------------------------- FASTQ parsing ----------------------------- #

def parse_fastq(filename: str | Path) -> List[Tuple[str, str]]:
    """Parse a FASTQ file and return (sequence, label) tuples.

    The "label" is the token immediately following '@' on the header line,
    up to the first whitespace. Order of records is preserved.

    Args:
        filename: Path to a FASTQ file (.fastq).

    Returns:
        A list of (sequence, label) tuples in file order.

    Raises:
        ValueError: If the FASTQ is malformed (e.g., bad header or incomplete record).
        OSError: If the file cannot be opened.
    """
    path = Path(filename)
    records: List[Tuple[str, str]] = []

    with path.open("r", encoding="utf-8") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            if not header.startswith("@"):
                raise ValueError("Malformed FASTQ: header line does not start with '@'.")

            label = header[1:].strip().split()[0]

            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not (seq and plus and qual):
                raise ValueError("Malformed FASTQ: incomplete record.")

            # strip trailing newlines/spaces
            seq = seq.strip()
            plus = plus.strip()
            qual = qual.strip()

            if not plus.startswith("+"):
                raise ValueError("Malformed FASTQ: third line must start with '+'.")

            records.append((seq, label))
    return records


def labels_longest_sequences(seq_label_list: List[Tuple[str, str]]) -> List[str]:
    """Return labels of sequences with maximum length.

    Preserves the input order for ties.

    Args:
        seq_label_list: List of (sequence, label) tuples.

    Returns:
        A list of labels corresponding to all sequences that share
        the maximum length. Returns an empty list if input is empty.
    """
    if not seq_label_list:
        return []
    max_len = max(len(seq) for seq, _ in seq_label_list)
    return [label for seq, label in seq_label_list if len(seq) == max_len]


# ----------------------------- ZIP handling ----------------------------- #

def unzip_to_script_dir(zip_path: Path) -> Path:
    """Extract a ZIP to the script directory and return the .fastq path.

    On Linux/macOS, this uses the external `unzip` command (quietly, no output).
    On Windows, it uses Python's `zipfile` module. The extracted file is expected
    to have the same base name as the zip, with a `.fastq` extension.

    Args:
        zip_path: Path to the input ZIP archive.

    Returns:
        Path to the extracted `.fastq` file located in the script directory.

    Raises:
        SystemExit: If the OS is unsupported or the file is invalid.
        CalledProcessError: If the external `unzip` command fails (Linux/macOS).
        FileNotFoundError: If the expected `.fastq` is not found after extraction.
    """
    script_dir = Path(__file__).resolve().parent
    expected_fastq = script_dir / f"{zip_path.stem}.fastq"

    system = platform.system()

    if system in {"Linux", "Darwin"}:
        try:
            subprocess.run(
                ["unzip", "-qq", "-o", str(zip_path), "-d", str(script_dir)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError:
            sys.exit("Error: invalid file")
    elif system == "Windows":
        try:
            with zipfile.ZipFile(zip_path, "r") as zf:
                zf.extractall(script_dir)
        except zipfile.BadZipFile:
            sys.exit("Error: invalid file")
    else:
        sys.exit(f"Error: unsupported operating system '{system}'")

    if not expected_fastq.is_file():
        raise FileNotFoundError(f"Expected extracted FASTQ not found ({expected_fastq.name})")

    return expected_fastq


# ----------------------------- Main logic ----------------------------- #

def main() -> None:
    """CLI entry point: handle input, extract if needed, and print results.

    Exits with:
        - Usage line if the number of arguments is incorrect.
        - "Error: invalid file" if the provided path is not a file.

    Side Effects:
        Prints the labels of the longest sequence(s) to stdout.
    """
    if len(sys.argv) != 2:
        sys.exit(f"Usage: python3 {Path(sys.argv[0]).name} <FILENAME.fastq or FILENAME.zip>")

    input_arg = Path(sys.argv[1])

    if not input_arg.is_file():
        sys.exit("Error: invalid file")

    if input_arg.suffix.lower() == ".zip":
        fastq_path = unzip_to_script_dir(input_arg)
    else:
        fastq_path = input_arg

    try:
        records = parse_fastq(fastq_path)
    except (OSError, ValueError) as e:
        sys.exit(f"Error: {e}")

    for label in labels_longest_sequences(records):
        print(label)


if __name__ == "__main__":
    main()
