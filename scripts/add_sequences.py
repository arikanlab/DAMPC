import argparse
import os
from typing import Optional
import pandas as pd
from Bio import SeqIO
from itertools import islice

def normalize_id(x: str, splitter: Optional[str]) -> str:
    if splitter is None:
        return x
    return x.split(splitter, 1)[0]

def main():
    parser = argparse.ArgumentParser(
        description="Join sequences from a FASTA file into a (tabular) input by matching IDs, using low memory."
    )
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("input_file", help="Path to the input TXT/CSV (must contain the ID column)")
    parser.add_argument("output_file", help="Path to the output file")
    parser.add_argument("--id-column", default="ID", help="Input ID column name (default: ID)")
    parser.add_argument("--seq-column", default="Sequence", help="Output sequence column name (default: Sequence)")
    parser.add_argument("--sep", default="\t", help=r"Field separator for input/output (default: '\t')")
    parser.add_argument("--chunksize", type=int, default=100_000, help="Rows per chunk (default: 100000)")
    parser.add_argument("--fasta-id-split", default=" ",
                        help="Use only the part of the FASTA ID before this string (default: space). "
                             "Use '' to disable and take the full FASTA ID.")
    parser.add_argument("--encoding", default="utf-8", help="Encoding for input/output tables (default: utf-8)")
    parser.add_argument("--verbose", action="store_true", help="Print extra diagnostics")
    args = parser.parse_args()

    splitter = None if args.fasta_id_split == "" else args.fasta_id_split

    # 0) Basic input checks
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Input file not found: {args.input_file}")
    if not os.path.exists(args.fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {args.fasta_file}")

    # 1) Build on-disk index for FASTA (low-memory)
    if args.verbose:
        print(f"[INFO] Indexing FASTA: {args.fasta_file}")
    fasta_index = SeqIO.index(args.fasta_file, "fasta")  # dict-like, lazy

    # Optionally show a few example FASTA IDs
    if args.verbose:
        sample_ids = list(islice(fasta_index.keys(), 5))
        print(f"[INFO] Example FASTA IDs (first 5): {sample_ids}")

    # 2) Stream input, map sequences, and write output incrementally
    if os.path.exists(args.output_file):
        os.remove(args.output_file)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_file)) or ".", exist_ok=True)

    header_written = False
    total_rows = 0
    matched_rows = 0

    # First chunk just to validate columns early (and to get dtype=str consistently)
    reader = pd.read_csv(
        args.input_file, sep=args.sep, chunksize=args.chunksize, dtype=str, encoding=args.encoding
    )

    for i, chunk in enumerate(reader, start=1):
        if args.id_column not in chunk.columns:
            raise ValueError(f"Input file must contain '{args.id_column}' column. Columns found: {list(chunk.columns)}")

        # Normalize IDs (both input and FASTA side)
        # We map by applying normalize_id on input IDs, and looking up in fasta_index by normalized key.
        norm_ids = chunk[args.id_column].astype(str).fillna("").map(lambda x: normalize_id(x, splitter))

        # Create sequence column
        # Important: fasta_index keys are based on record.id (not header line). If your IDs contain spaces etc.,
        # normalize_id() makes them comparable.
        seqs = []
        mcount = 0
        for rid in norm_ids:
            if rid and rid in fasta_index:
                mcount += 1
                seqs.append(str(fasta_index[rid].seq))
            else:
                seqs.append(None)
        matched_rows += mcount

        # Insert sequence column at position 1 (right after ID)
        chunk.insert(1, args.seq_column, pd.Series(seqs, index=chunk.index))

        # Write out
        chunk.to_csv(
            args.output_file,
            sep=args.sep,
            index=False,
            mode="a",
            header=not header_written,
            encoding=args.encoding,
            lineterminator="\n",
        )
        header_written = True
        total_rows += len(chunk)

        if args.verbose:
            print(f"[INFO] Wrote chunk {i}: rows={len(chunk)}, matched={mcount}")

    if args.verbose:
        print(f"[INFO] Done. Total rows: {total_rows}, matched rows: {matched_rows}")
    print(f"✅ Output saved to: {args.output_file}")

if __name__ == "__main__":
    main()
