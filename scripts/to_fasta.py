#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
from typing import Iterator, List, Tuple, Optional, Union

CANDIDATE_DELIMS = [",", "\t", ";", "|"]

def parse_col_arg(val: str) -> Union[int, str]:
    """Return int if numeric, else lowercase string."""
    s = val.strip()
    return int(s) if s.isdigit() else s.lower()

def sniff_delimiter(path: str, user_delim: Optional[str]) -> Optional[str]:
    """
    Determine delimiter.
    Priority:
      1) user-provided --delimiter
      2) by file extension (.csv -> ',', .tsv -> '\\t')
      3) csv.Sniffer over a sample
      4) heuristic count of candidate delimiters
      5) None => whitespace split
    """
    if user_delim:
        return user_delim

    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        return ","
    if ext == ".tsv":
        return "\t"

    # Read a small sample for detection (doesn't assume where header is)
    sample_lines: List[str] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for _ in range(100):
            line = f.readline()
            if not line:
                break
            if line.strip():
                sample_lines.append(line)

    sample_text = "".join(sample_lines)

    # Try csv.Sniffer first
    if sample_text:
        try:
            dialect = csv.Sniffer().sniff(sample_text, delimiters="".join(CANDIDATE_DELIMS))
            return dialect.delimiter
        except Exception:
            pass

    # Fallback: count candidate delimiters
    counts = {d: 0 for d in CANDIDATE_DELIMS}
    for line in sample_lines:
        for d in CANDIDATE_DELIMS:
            counts[d] += line.count(d)
    if counts and max(counts.values()) > 0:
        return max(counts, key=counts.get)

    # Final fallback: whitespace
    return None

def row_reader(path: str, delim: Optional[str]) -> Iterator[List[str]]:
    """
    Yield rows as list of fields.
    If delim is None, split on whitespace.
    Keeps blank-line skipping minimal; returns non-empty rows only.
    """
    if delim is None:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                yield re.split(r"\s+", line)
        return

    with open(path, "r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        for row in reader:
            if not row:
                continue
            yield [cell.strip() for cell in row]

def resolve_indices(
    header: List[str],
    id_col: Union[int, str],
    seq_col: Union[int, str]
) -> Tuple[int, int]:
    """
    Map provided column specifiers to indices using the header if names given.
    """
    if isinstance(id_col, int) and isinstance(seq_col, int):
        return id_col, seq_col

    lower = [h.strip().lower() for h in header]

    def idx_for(col):
        if isinstance(col, int):
            return col
        try:
            return lower.index(col)  # col is lowercase name
        except ValueError:
            raise SystemExit(f"Column '{col}' not found in header: {lower}")

    return idx_for(id_col), idx_for(seq_col)

def extract_records(
    path: str,
    id_col: Union[int, str],
    seq_col: Union[int, str],
    delimiter_opt: Optional[str],
    uppercase: bool,
    header_row: int,
) -> Iterator[Tuple[str, str]]:
    """
    Extract (id, seq) pairs.
    header_row: 1-based row index of header; 0 means 'no header'.
    """
    if header_row < 0:
        raise SystemExit("--header-row must be >= 0")

    delim = sniff_delimiter(path, delimiter_opt)
    reader = row_reader(path, delim)

    current_row_num = 0
    header: Optional[List[str]] = None

    # Advance to header row (if any)
    for row in reader:
        current_row_num += 1
        if header_row == 0:
            # No header mode: first non-empty row is data
            header = None
            # But we need indices: require integer specs
            if not (isinstance(id_col, int) and isinstance(seq_col, int)):
                raise SystemExit("With --header-row 0 (no header), you must use numeric indices for --id-col and --seq-col.")
            # Process this row and then break to main loop with known indices
            first_data_row = row
            id_idx, seq_idx = id_col, seq_col
            # Yield first row then continue with remaining rows
            rid = str(first_data_row[id_idx]).strip().strip('"').strip("'") if len(first_data_row) > id_idx else ""
            seq = str(first_data_row[seq_idx]).strip().strip('"').strip("'") if len(first_data_row) > seq_idx else ""
            if uppercase:
                seq = seq.upper()
            if rid and seq:
                yield rid, seq
            break
        else:
            if current_row_num == header_row:
                header = row
                break
            else:
                continue
    else:
        # File exhausted before reaching header/data
        raise SystemExit(f"No data found in {path}")

    # If we’re here and header_row > 0, we have a header
    id_idx, seq_idx = resolve_indices(header, id_col, seq_col)

    # Remaining rows are data
    for row in reader:
        if len(row) <= max(id_idx, seq_idx):
            continue
        rid = str(row[id_idx]).strip().strip('"').strip("'")
        seq = str(row[seq_idx]).strip().strip('"').strip("'")
        if uppercase:
            seq = seq.upper()
        if rid and seq:
            yield rid, seq

def wrap_seq(seq: str, width: int) -> str:
    if width <= 0:
        return seq + "\n"
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width)) + "\n"

def write_fasta(pairs: Iterator[Tuple[str, str]], out_path: str, wrap: int) -> int:
    count = 0
    with open(out_path, "w", encoding="utf-8") as out:
        for rid, seq in pairs:
            out.write(f">{rid}\n")
            out.write(wrap_seq(seq, wrap))
            count += 1
    return count

def main():
    ap = argparse.ArgumentParser(
        description="Convert CSV/TSV/TXT into FASTA with user-specified ID/sequence columns. Supports header on any line."
    )
    ap.add_argument("input", help="Input file (.csv, .tsv, .txt)")
    ap.add_argument("-o", "--out", required=True, help="Output FASTA path")
    ap.add_argument("--id-col", required=True, help="ID column (name or 0-based index)")
    ap.add_argument("--seq-col", required=True, help="Sequence column (name or 0-based index)")
    ap.add_argument("--header-row", type=int, default=1,
                    help="1-based header row index; use 0 if there is no header. Default: 1")
    ap.add_argument("--delimiter", help="Override delimiter (e.g., ',' or '\\t'). If omitted, auto-detect.")
    ap.add_argument("--uppercase", action="store_true", help="Uppercase sequences")
    ap.add_argument("--wrap", type=int, default=60, help="FASTA line wrap width (0 = no wrap). Default: 60")
    args = ap.parse_args()

    id_col = parse_col_arg(args.id_col)
    seq_col = parse_col_arg(args.seq_col)

    pairs = list(extract_records(
        args.input, id_col, seq_col,
        delimiter_opt=args.delimiter,
        uppercase=args.uppercase,
        header_row=args.header_row
    ))
    if not pairs:
        sys.exit("No valid (id, sequence) rows found. Check column mapping / header row.")

    n = write_fasta(pairs, args.out, args.wrap)
    print(f"Wrote {n} records to {args.out}")

# reuse parse_col_arg in main
def parse_col_arg(val: str) -> Union[int, str]:
    s = val.strip()
    return int(s) if s.isdigit() else s.lower()

if __name__ == "__main__":
    main()
