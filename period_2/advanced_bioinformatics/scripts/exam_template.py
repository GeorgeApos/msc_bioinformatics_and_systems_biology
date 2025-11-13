#!/usr/bin/env python3
"""
Abstract Exam Template (no 3rd-party libs)

What you get:
- Basic CLI with common flags (input, ref, output, threads, dry-run, force).
- Ready-made external command runner with caching to <output>.
- Minimal parsing skeleton (line-by-line; header-aware).
- Compute + output section (print or write TSV).

Typical usage:
  python exam.py -i input.fa -r ref.fa -o out.txt
  python exam.py -i reads.fq -o out.tsv --dry-run           # show the command only
  python exam.py -i data.txt -o parsed.tsv --force          # overwrite cached output
"""

from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Dict, Iterable, Optional
import argparse
import subprocess
import sys

# -------------------------- utils -------------------------- #

def die(msg: str, code: int = 1) -> None:
    print(msg, file=sys.stderr)
    sys.exit(code)

def ensure_parent(path: Path) -> None:
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)

# ---------------------- argument parsing ------------------- #

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Abstract exam runner + parser scaffold."
    )
    p.add_argument("-i", "--input", required=True, help="Primary input file.")
    p.add_argument("-r", "--ref", help="Optional reference/support file.")
    p.add_argument("-o", "--output", required=True, help="Output file to write.")
    p.add_argument("-t", "--threads", type=int, default=1, help="Threads/CPUs to use.")
    p.add_argument("--dry-run", action="store_true", help="Print the command; do not execute.")
    p.add_argument("--force", action="store_true", help="Re-run even if output exists.")
    p.add_argument("--quiet", action="store_true", help="Less stderr chatter.")
    return p

# ------------------- command construction ------------------ #

def build_command(args: argparse.Namespace) -> List[str]:
    """
    TODO: Replace with the real command they give you in the exam.
    Keep each token a separate list element (no shell=True).
    Example placeholder below.
    """
    cmd = [
        "toolname",                 # <- change to real executable
        "--input", args.input,      # or args.ref, etc.
        "--threads", str(args.threads),
        # ... add flags as required in the exam ...
    ]
    return cmd

# ------------------- command execution --------------------- #

def run_command_once(cmd: List[str], outfile: Path, dry_run: bool, force: bool, quiet: bool) -> None:
    """
    Runs the external command once to produce 'outfile'.
    Skips re-run if 'outfile' exists and is non-empty, unless --force.
    """
    if outfile.exists() and outfile.stat().st_size > 0 and not force:
        if not quiet:
            print(f"[cache] Using existing {outfile}", file=sys.stderr)
        return

    if dry_run:
        print("[dry-run] Command:", " ".join(cmd), file=sys.stderr)
        return

    ensure_parent(outfile)

    if not quiet:
        print("[run] ", " ".join(cmd), file=sys.stderr)

    try:
        with outfile.open("w", encoding="utf-8") as out:
            # If the real tool writes to files by itself, you may prefer capture_output=True here.
            res = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, check=False)
    except FileNotFoundError:
        die(f"Error: command not found: {cmd[0]}")
    if res.returncode != 0:
        err = (res.stderr or "").strip()
        die(f"External tool failed (exit {res.returncode}).\n{err if err else '(no stderr)'}")

# --------------------- parsing scaffold -------------------- #

@dataclass
class ParsedRow:
    # Fill/rename these fields for the exam
    col1: str
    col2: int
    col3: float

def parse_output_file(path: Path) -> List[ParsedRow]:
    """
    Minimal, robust parser:
    - Skips empty/comment lines.
    - Detects/ignores a header if first token is non-numeric where numeric is expected.
    - Splits on tab by default (adjust as needed).
    """
    rows: List[ParsedRow] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")   # adjust delimiter if needed
            # ---- TODO: map cols → ParsedRow fields properly ----
            # Example mapping with basic guards:
            try:
                # If first line is a header, this will fail and we continue.
                c1 = cols[0]
                c2 = int(cols[1])
                c3 = float(cols[2])
            except (IndexError, ValueError):
                # Probably a header or malformed line; skip gracefully.
                continue
            rows.append(ParsedRow(c1, c2, c3))
    return rows

# ------------------ compute / transform step ---------------- #

@dataclass
class Result:
    n_rows: int
    sum_col2: int
    mean_col3: float
    top_items: List[Tuple[str, int]]  # e.g., top by col2

def compute_results(parsed: List[ParsedRow], top_k: int = 5) -> Result:
    if not parsed:
        return Result(0, 0, 0.0, [])

    n = len(parsed)
    sum2 = sum(p.col2 for p in parsed)
    mean3 = sum(p.col3 for p in parsed) / n

    # Simple “top K by col2”
    items = sorted(((p.col1, p.col2) for p in parsed), key=lambda x: x[1], reverse=True)[:top_k]

    return Result(n_rows=n, sum_col2=sum2, mean_col3=mean3, top_items=items)

# ------------------------- output section ------------------- #

def print_report(res: Result) -> None:
    """
    Human-readable summary (stdout). Keep it deterministic and concise.
    """
    print("# Summary")
    print(f"rows\t{res.n_rows}")
    print(f"sum_col2\t{res.sum_col2}")
    print(f"mean_col3\t{res.mean_col3:.3f}")
    print()
    print("# Top items by col2")
    print("rank\titem\tcol2")
    for i, (k, v) in enumerate(res.top_items, start=1):
        print(f"{i}\t{k}\t{v}")

def write_tsv(rows: List[ParsedRow], out: Path) -> None:
    """
    Optional: write normalized rows as TSV for grading scripts.
    """
    ensure_parent(out)
    with out.open("w", encoding="utf-8") as fh:
        fh.write("col1\tcol2\tcol3\n")
        for r in rows:
            fh.write(f"{r.col1}\t{r.col2}\t{r.col3}\n")

# ----------------------------- main ------------------------ #

def main() -> None:
    args = build_parser().parse_args()
    inpath = Path(args.input).resolve()
    refpath = Path(args.ref).resolve() if args.ref else None
    outpath = Path(args.output).resolve()

    # 1) Build the exact command they specify in the exam.
    cmd = build_command(args)

    # 2) Run once (with caching), or just show it in --dry-run mode.
    run_command_once(cmd, outpath, dry_run=args.dry_run, force=args.force, quiet=args.quiet)

    # If we were in dry-run and no file was produced, stop here.
    if args.dry_run:
        return

    # 3) Parse the produced output (or you can parse another file if the tool writes elsewhere).
    parsed = parse_output_file(outpath)

    # 4) Compute simple results (edit this quickly to match the question).
    res = compute_results(parsed, top_k=5)

    # 5) Print a neat report (stdout). Optionally also write a normalized TSV.
    print_report(res)
    # write_tsv(parsed, Path("normalized.tsv"))  # uncomment if needed

if __name__ == "__main__":
    main()