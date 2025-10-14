import argparse
import os
from typing import Optional, Dict, Set, Iterable
import pandas as pd
from Bio import SeqIO

# -------- Normalization helpers --------
def norm_seq(s: str, upper: bool, strip_ws: bool) -> str:
    if s is None:
        return ""
    s2 = str(s)
    if strip_ws:
        s2 = "".join(s2.split())  # remove all whitespace
    if upper:
        s2 = s2.upper()
    return s2

def norm_fasta_id(rid: str, splitter: Optional[str]) -> str:
    if splitter is None:
        return rid
    return rid.split(splitter, 1)[0]

# -------- Collect needed sequences from CSV (streaming) --------
def collect_needed_sequences(input_path: str, seq_col: str, sep: str, encoding: str,
                             chunksize: int, upper: bool, strip_ws: bool) -> Set[str]:
    needed: Set[str] = set()
    reader = pd.read_csv(input_path, sep=sep, chunksize=chunksize, dtype=str, encoding=encoding)
    for chunk in reader:
        if seq_col not in chunk.columns:
            raise ValueError(f"Input file must contain '{seq_col}' column. Columns: {list(chunk.columns)}")
        seqs = (norm_seq(x, upper, strip_ws) for x in chunk[seq_col].dropna().tolist())
        needed.update(s for s in seqs if s)
    return needed

# -------- Build map Sequence -> ID by scanning FASTA once --------
def build_seq_to_id_map(fasta_path: str, needed: Set[str], fasta_id_split: Optional[str],
                        upper: bool, strip_ws: bool) -> Dict[str, str]:
    seq2id: Dict[str, str] = {}
    if not needed:
        return seq2id
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq_norm = norm_seq(str(rec.seq), upper, strip_ws)
        if seq_norm in needed and seq_norm not in seq2id:
            seq2id[seq_norm] = norm_fasta_id(rec.id, fasta_id_split)
    return seq2id

# -------- Stream CSV again, write IDs (or leave blank when missing) --------
def write_with_ids(input_path: str, output_path: str, id_col: str, seq_col: str, sep: str, encoding: str,
                   chunksize: int, seq2id: Dict[str, str], upper: bool, strip_ws: bool, insert_pos: int) -> None:
    if os.path.exists(output_path):
        os.remove(output_path)
    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)

    header_written = False
    reader = pd.read_csv(input_path, sep=sep, chunksize=chunksize, dtype=str, encoding=encoding)
    for chunk in reader:
        if seq_col not in chunk.columns:
            raise ValueError(f"Input file must contain '{seq_col}' column. Columns: {list(chunk.columns)}")
        # map each sequence to ID (missing -> None)
        ids = []
        for s in chunk[seq_col].astype(str):
            key = norm_seq(s, upper, strip_ws)
            ids.append(seq2id.get(key, None))
        chunk.insert(max(0, insert_pos), id_col, ids)
        chunk.to_csv(output_path, sep=sep, index=False, mode="a", header=not header_written,
                     encoding=encoding, lineterminator="\n")
        header_written = True

def main():
    p = argparse.ArgumentParser(
        description="Match sequences in a results file to a FASTA and write IDs for matches; skip (leave blank) when not found."
    )
    p.add_argument("fasta_file", help="Path to the FASTA file")
    p.add_argument("input_file", help="Path to the input results file (CSV/TSV)")
    p.add_argument("output_file", help="Path to the output file")
    p.add_argument("--seq-column", default="Sequence", help="Column name in input that holds sequences (default: Sequence)")
    p.add_argument("--id-column", default="ID", help="Name of the ID column to add (default: ID)")
    p.add_argument("--insert-pos", type=int, default=0, help="Insert position for the ID column (default: 0 = as first column)")
    p.add_argument("--sep", default=",", help=r"Field separator for input/output (default: ','). Use '\t' for TSV.")
    p.add_argument("--encoding", default="utf-8", help="Encoding for input/output (default: utf-8)")
    p.add_argument("--chunksize", type=int, default=200_000, help="Rows per processing chunk (default: 200000)")
    p.add_argument("--upper", action="store_true", help="Upper-case both CSV and FASTA sequences before matching")
    p.add_argument("--strip-ws", action="store_true", help="Remove all whitespace (spaces, tabs, newlines) before matching")
    p.add_argument("--fasta-id-split", default=" ",
                   help="Use only the part of FASTA ID before this string (default: space). Use '' to disable.")
    p.add_argument("--verbose", action="store_true", help="Print progress information")
    args = p.parse_args()

    splitter = None if args.fasta_id_split == "" else args.fasta_id_split

    if args.verbose:
        print("[INFO] Collecting unique sequences from input (streaming)...")
    needed = collect_needed_sequences(args.input_file, args.seq_column, args.sep, args.encoding,
                                      args.chunksize, args.upper, args.strip_ws)
    if args.verbose:
        print(f"[INFO] Unique sequences needed: {len(needed):,}")

    if args.verbose:
        print("[INFO] Scanning FASTA to build Sequence->ID map...")
    seq2id = build_seq_to_id_map(args.fasta_file, needed, splitter, args.upper, args.strip_ws)
    if args.verbose:
        print(f"[INFO] Matches found in FASTA: {len(seq2id):,}")

    if args.verbose and len(seq2id) < len(needed):
        print(f"[WARN] {len(needed) - len(seq2id):,} sequences from input were not found in FASTA (they will be left blank).")

    if args.verbose:
        print("[INFO] Writing output (streaming)...")
    write_with_ids(args.input_file, args.output_file, args.id_column, args.seq_column, args.sep,
                   args.encoding, args.chunksize, seq2id, args.upper, args.strip_ws, args.insert_pos)

    print(f"✅ Output saved to: {args.output_file}")

if __name__ == "__main__":
    main()
