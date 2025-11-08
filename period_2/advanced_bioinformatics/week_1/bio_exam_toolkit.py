#!/usr/bin/env python3
"""
Bio Exam Toolkit — minimal helpers for parsing & basic bioinformatics tasks.

What’s inside (functions):
- Generic
  * length_stats, tab_print
- FASTA
  * parse_fasta, write_fasta, wrap_fasta
- FASTQ
  * parse_fastq, detect_quality_offset, convert_quality_line, average_q_per_pos
- GenBank (flatfile, no external libs)
  * parse_genbank
- DNA/RNA & proteins
  * gc_content, reverse_complement, dna_to_rna, translate_rna, find_orfs_6frames
- K-mers
  * extract_kmers, kmer_stats
- SMILES
  * smiles_ring_count
- Tools
  * run_cmd, run_fastq_quality_trimmer
"""

from __future__ import annotations
from pathlib import Path
from typing import Iterable, Iterator, List, Tuple, Dict
import subprocess
import sys
import os
import csv

# --------------------- Generic utilities ---------------------

def length_stats(seqs: Iterable[str]) -> Tuple[int, int, float]:
    """Return (min_len, max_len, avg_len) for an iterable of sequences."""
    lens = [len(s) for s in seqs]
    if not lens:
        raise ValueError("No sequences provided.")
    return min(lens), max(lens), sum(lens)/len(lens)

def tab_print(*cols) -> None:
    """Print columns separated by tabs (no trailing spaces)."""
    print("\t".join(str(c) for c in cols))

# --------------------- FASTA parsing & writing ---------------------

def parse_fasta(src: Iterable[str] | Path) -> Dict[str, str]:
    """Parse FASTA and return dict{id: sequence} (concatenates wrapped lines)."""
    if isinstance(src, (str, Path)):
        fh: Iterator[str]
        with open(src, "r", encoding="utf-8") as fh:
            return parse_fasta(fh)

    seqs: Dict[str, str] = {}
    label = None
    parts: List[str] = []
    for raw in src:
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if label is not None:
                seqs[label] = "".join(parts)
            label = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)
    if label is not None:
        seqs[label] = "".join(parts)
    return seqs

def wrap_fasta(seq: str, width: int = 60) -> List[str]:
    """Wrap a sequence into fixed-width lines."""
    return [seq[i:i+width] for i in range(0, len(seq), width)]

def write_fasta(records: Dict[str, str], out_fn: str | Path, width: int = 60) -> None:
    """Write dict{id: seq} to a FASTA file, wrapped at `width`."""
    with open(out_fn, "w", encoding="utf-8") as f:
        for rid, seq in records.items():
            f.write(f">{rid}\n")
            for line in wrap_fasta(seq, width):
                f.write(line + "\n")

# --------------------- FASTQ parsing & qualities ---------------------

def parse_fastq(path: str | Path) -> List[Tuple[str, str, str]]:
    """Return list of (label, sequence, quality) from a FASTQ file."""
    out: List[Tuple[str, str, str]] = []
    p = Path(path)
    with p.open("r", encoding="utf-8") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            if not header.startswith("@"):
                raise ValueError("Malformed FASTQ: header must start with '@'.")
            label = header[1:].strip().split()[0]

            seq  = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not (seq and plus and qual):
                raise ValueError("Malformed FASTQ: incomplete record.")
            if not plus.startswith("+"):
                raise ValueError("Malformed FASTQ: third line must start with '+'.")
            out.append((label, seq.strip(), qual.strip()))
    return out

def detect_quality_offset(qual_lines: List[str]) -> int:
    """Detect Phred offset (33 or 64) from raw quality lines."""
    if not qual_lines:
        return 33
    min_code = min(ord(ch) for q in qual_lines for ch in q)
    return 64 if min_code >= 64 else 33

def convert_quality_line(qstr: str, q_offset: int = 64) -> List[int]:
    """Convert ASCII quality string to list of ints; clamp to [0, 41]."""
    out: List[int] = []
    for ch in qstr:
        q = ord(ch) - q_offset
        if q < 0:   q = 0
        if q > 41:  q = 41
        out.append(q)
    return out

def average_q_per_pos(qualities: List[List[int]]) -> List[float]:
    """Average quality at each base position across reads (ragged-aware)."""
    if not qualities:
        return []
    max_len = max(len(q) for q in qualities)
    sums = [0]*max_len
    counts = [0]*max_len
    for q in qualities:
        for i, v in enumerate(q):
            sums[i] += v
            counts[i] += 1
    return [sums[i]/counts[i] if counts[i] else 0.0 for i in range(max_len)]

# --------------------- GenBank (flat) parsing ---------------------

def parse_genbank(path: str | Path) -> List[Tuple[str, str, str]]:
    """
    Parse GenBank flatfile into list of (accession, organism, sequence).
    - accession: first token on ACCESSION line
    - organism:  text on the '  ORGANISM' line
    - sequence:  letters from ORIGIN..'//', uppercase, letters only
    """
    records: List[Tuple[str, str, str]] = []
    accession = None
    organism  = None
    seq_chunks: List[str] = []
    in_rec = False
    in_seq = False

    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("LOCUS"):
                if in_rec and accession and organism:
                    seq = "".join(seq_chunks)
                    records.append((accession, organism, seq))
                accession = None
                organism = None
                seq_chunks = []
                in_rec = True
                in_seq = False
                continue

            if not in_rec:
                continue

            if line.startswith("ACCESSION"):
                parts = line.split()
                if len(parts) >= 2:
                    accession = parts[1]
                continue

            if line.startswith("  ORGANISM"):
                organism = line.split("ORGANISM", 1)[1].strip() if "ORGANISM" in line else ""
                continue

            if line.startswith("ORIGIN"):
                in_seq = True
                continue

            if line.strip() == "//":
                seq = "".join(seq_chunks)
                records.append((accession or "", organism or "", seq))
                accession = organism = None
                seq_chunks = []
                in_rec = in_seq = False
                continue

            if in_seq:
                letters = "".join(ch for ch in line if ch.isalpha())
                if letters:
                    seq_chunks.append(letters.upper())

    return records

# --------------------- DNA/RNA & protein helpers ---------------------

def gc_content(seq: str) -> float:
    """GC percentage (0..100). Counts only A/C/G/T in denominator."""
    if not seq:
        return 0.0
    s = seq.upper()
    atgc = sum(s.count(x) for x in "ACGT")
    if atgc == 0:
        return 0.0
    return (s.count("G") + s.count("C")) / atgc * 100.0

def reverse_complement(dna: str) -> str:
    """Reverse complement of DNA (A<->T, C<->G)."""
    comp = {"A":"T","T":"A","C":"G","G":"C"}
    dna = dna.upper()
    return "".join(comp[b] for b in reversed(dna))

def dna_to_rna(dna: str) -> str:
    """DNA -> RNA by T->U, uppercase."""
    return dna.upper().replace("T", "U")

def _genetic_code_rna() -> Dict[str, str]:
    """Standard genetic code as RNA codons -> AA (single-letter)."""
    first  = 16 * 'U' + 16 * 'C' + 16 * 'A' + 16 * 'G'
    second = (4 * 'U' + 4 * 'C' + 4 * 'A' + 4 * 'G') * 4
    third  = 16 * 'UCAG'
    prot   = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    return {first[i] + second[i] + third[i]: prot[i] for i in range(len(first))}
_CODE_RNA = _genetic_code_rna()

def translate_rna(rna: str) -> str:
    """Translate RNA from start to first stop (if present); ignores partial codons at tail."""
    aa = []
    for i in range(0, len(rna) - (len(rna) % 3), 3):
        sym = _CODE_RNA.get(rna[i:i+3])
        if sym is None:
            break
        if sym == "*":
            break
        aa.append(sym)
    return "".join(aa)

def find_orfs_6frames(dna: str) -> List[str]:
    """Find translations of all ORFs in 6 frames (ATG..TAA|TAG|TGA)."""
    dna = dna.upper()
    stops = {"TAA","TAG","TGA"}
    def frame_orfs(s: str, offset: int) -> List[str]:
        out: List[str] = []
        n = len(s)
        i = offset
        while i <= n-3:
            if s[i:i+3] == "ATG":
                j = i+3
                while j <= n-3:
                    cod = s[j:j+3]
                    if cod in stops:
                        out.append(translate_rna(dna_to_rna(s[i:j+3])))
                        break
                    j += 3
            i += 3
        return out

    rc = reverse_complement(dna)
    proteins: List[str] = []
    for f in (0,1,2):
        proteins.extend(frame_orfs(dna, f))
    for f in (0,1,2):
        proteins.extend(frame_orfs(rc, f))
    return proteins

# --------------------- K-mers ---------------------

def extract_kmers(seqs: Dict[str, str], k: int) -> Dict[str, int]:
    """Return {kmer: count} across all sequences (no RC collapsing)."""
    counts: Dict[str, int] = {}
    for seq in seqs.values():
        n = len(seq)
        if n < k:
            continue
        for i in range(n - k + 1):
            kmer = seq[i:i+k]
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts

def kmer_stats(kmer_table: Dict[str, int]) -> Tuple[int,int,int,int]:
    """Return (unique_count, distinct, total, max_count)."""
    distinct = len(kmer_table)
    total = sum(kmer_table.values())
    unique_count = sum(1 for c in kmer_table.values() if c == 1)
    max_count = max(kmer_table.values()) if kmer_table else 0
    return unique_count, distinct, total, max_count

# --------------------- SMILES rings ---------------------

def smiles_ring_count(smiles: str) -> int:
    """Count rings in a SMILES by counting un-bracketed digits and halving."""
    in_bracket = False
    digits = 0
    for ch in smiles:
        if ch == "[":
            in_bracket = True
        elif ch == "]":
            in_bracket = False
        elif not in_bracket and ch.isdigit():
            digits += 1
    return digits // 2

# --------------------- External tools ---------------------

def run_cmd(cmd: List[str]) -> subprocess.CompletedProcess:
    """Run a command and return CompletedProcess (raises on non-zero)."""
    return subprocess.run(cmd, check=True, capture_output=True, text=True)

def run_fastq_quality_trimmer(input_fq: Path, output_fq: Path, threshold: int, q_offset: int) -> None:
    """Wrap FASTX 'fastq_quality_trimmer' (-t, -Q, -i, -o)."""
    if not shutil_which("fastq_quality_trimmer"):
        sys.exit("Error: fastq_quality_trimmer not found in PATH")
    cmd = ["fastq_quality_trimmer", "-t", str(threshold), "-Q", str(q_offset),
           "-i", str(input_fq), "-o", str(output_fq)]
    try:
        res = subprocess.run(cmd, check=False, capture_output=True, text=True)
    except Exception as e:
        sys.exit(f"Error running fastq_quality_trimmer: {e}")
    if res.returncode != 0:
        msg = res.stderr.strip() or "Unknown error from fastq_quality_trimmer"
        sys.exit(f"Error running fastq_quality_trimmer: {msg}")
    if (not output_fq.is_file()) or os.path.getsize(output_fq) == 0:
        sys.exit("Error: trimmed output not created or empty")

def shutil_which(name: str) -> str | None:
    """Tiny replacement for shutil.which to avoid import if desired."""
    for p in os.environ.get("PATH", "").split(os.pathsep):
        cand = Path(p) / name
        if cand.exists() and os.access(cand, os.X_OK):
            return str(cand)
    return None

# --------------------- Generic utilities ---------------------

def safe_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default

def fmt3(x):  # three decimals
    return f"{x:.3f}"

def top_n(items, n, key=lambda x: x, reverse=True):
    return sorted(items, key=key, reverse=reverse)[:n]

# --------------------- FASTA parsing & writing ---------------------


def fasta_to_ordered(path: Path):
    """Return (ids, seqs_dict) preserving order; handy to keep chr order stable."""
    ids = []
    seqs = {}
    with open(path, "r", encoding="utf-8") as fh:
        cur = None
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                ids.append(cur)
                seqs[cur] = []
            else:
                seqs[cur].append(line.strip())
    for k in ids:
        seqs[k] = "".join(seqs[k])
    return ids, seqs

# --------------------- TSV / expression metrics ---------------------

def read_tsv(path, required_cols, delimiter="\t"):
    """Return (header_map, rows_as_dicts). Header match is case-sensitive."""
    with open(path, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split(delimiter)
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in required_cols if c not in idx]
        if missing:
            raise ValueError(f"Missing columns: {missing}")
        rows = []
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split(delimiter)
            row = {name: parts[i] for name, i in idx.items() if i < len(parts)}
            rows.append(row)
    return idx, rows

def gene_of_transcript(tx_id: str) -> str:
    """Return gene id by stripping at the first dot: AT1G00010.1 -> AT1G00010."""
    return tx_id.split(".", 1)[0]

def rpkm_from_counts_efflen(est_counts: float, eff_length: float, total_counts: float) -> float:
    """RPKM = 1e9 * counts / (eff_length * total_counts)."""
    if eff_length <= 0 or total_counts <= 0:
        return 0.0
    return 1e9 * est_counts / (eff_length * total_counts)

def tpm_from_rpkm(rpkm_vals):
    """Return list of TPM given list of RPKM values."""
    denom = sum(rpkm_vals)
    if denom == 0:
        return [0.0] * len(rpkm_vals)
    return [1e6 * r / denom for r in rpkm_vals]

def group_sum(items, key_fn, val_fn):
    """Generic sum-by-key: items -> {key: sum(values)}."""
    out = {}
    for it in items:
        k = key_fn(it)
        out[k] = out.get(k, 0.0) + float(val_fn(it))
    return out

# --------------------- External tools & intervals ---------------------

def run_tool_once(cmd: list[str], out_path: Path | None = None) -> None:
    """
    Run external command unless `out_path` exists & is non-empty.
    Print stderr on failure and exit(1).
    """
    if out_path and out_path.is_file() and os.path.getsize(out_path) > 0:
        return
    try:
        res = subprocess.run(cmd, check=False, capture_output=True, text=True)
    except Exception as e:
        print(f"Error starting tool: {e}", file=sys.stderr)
        sys.exit(1)
    if res.returncode != 0:
        print(res.stderr.strip() or "Unknown tool error", file=sys.stderr)
        sys.exit(1)

def run_lastz_general(query_fa: Path, ref_fa: Path, out_txt: Path):
    cmd = ["lastz", str(query_fa), str(ref_fa), "--format=general", "--output="+str(out_txt)]
    run_tool_once(cmd, out_txt)
    return out_txt

def parse_lastz_general(path: Path):
    """
    Return list of alignments as dicts with keys:
    'refName','refStart','refEnd' (1-based inclusive as stored in file).
    This parser expects a header line containing these fields.
    """
    with open(path, "r", encoding="utf-8") as fh:
        header = fh.readline().strip().split()
        idx = {name: i for i, name in enumerate(header)}
        need = ["refName","refStart","refEnd"]
        for n in need:
            if n not in idx:
                raise ValueError(f"LASTZ general header missing '{n}'")
        recs = []
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split()
            recs.append({
                "refName": parts[idx["refName"]],
                "refStart": int(parts[idx["refStart"]]),
                "refEnd": int(parts[idx["refEnd"]]),
            })
    return recs

def merge_intervals(intervals):
    """Merge [start, end) half-open intervals; input must be sorted by start."""
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(a, b) for a, b in merged]

def complement_intervals(covered, length):
    """
    Given merged covered intervals [start,end) on [0,length), return uncovered as list of [s,e).
    """
    res = []
    prev = 0
    for s, e in covered:
        if prev < s:
            res.append((prev, s))
        prev = max(prev, e)
    if prev < length:
        res.append((prev, length))
    return res

def run_kallisto_quant(index_fn: Path, reads_fn: Path, out_dir: Path,
                       frag_len=200, frag_sd=20):
    out_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = out_dir / "abundance.tsv"
    cmd = ["kallisto", "quant",
           "-i", str(index_fn),
           "--single", "-l", str(frag_len), "-s", str(frag_sd),
           "-o", str(out_dir),
           str(reads_fn)]
    run_tool_once(cmd, out_tsv)
    return out_tsv

def accept_args():
    """Accept command-line arguments.

    Returns:
        filenames: list of str, input filenames
    """

    if len(sys.argv) != 3:
        print("Usage: e.g. python3 my_exam.py velvet_15.fa chr3.fa")
        sys.exit(1)

    filenames = sys.argv[1:]
    return filenames