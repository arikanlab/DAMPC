import argparse

# Set of 20 standard amino acids
NATURAL_AAS = set("ARNDCEQGHILKMFPSTWYV")

def analyze_fasta(filepath):
    """
    Parse a FASTA file and return:
      - total number of sequences (including invalid)
      - dict mapping valid sequence -> first seen header
    A valid sequence contains only the 20 standard amino acids.
    """
    total = 0
    seq_to_header = {}
    current_seq = []
    current_header = None

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):  # header line
                if current_seq:
                    seq_str = "".join(current_seq).upper()
                    total += 1
                    if seq_str and all(ch in NATURAL_AAS for ch in seq_str):
                        seq_to_header.setdefault(seq_str, current_header)
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        # add the last sequence
        if current_seq:
            seq_str = "".join(current_seq).upper()
            total += 1
            if seq_str and all(ch in NATURAL_AAS for ch in seq_str):
                seq_to_header.setdefault(seq_str, current_header)

    return total, seq_to_header


def main():
    parser = argparse.ArgumentParser(
        description="Compare two FASTA files and extract common and unique sequences."
    )
    parser.add_argument("fasta1", help="Path to the first FASTA file")
    parser.add_argument("fasta2", help="Path to the second FASTA file")
    args = parser.parse_args()

    # Analysis
    total1, seqs1 = analyze_fasta(args.fasta1)
    total2, seqs2 = analyze_fasta(args.fasta2)

    set1 = set(seqs1.keys())
    set2 = set(seqs2.keys())

    # ---- Comparison ----
    common = set1 & set2
    only_second = set2 - set1

    print("\n=== Comparison ===")
    print(f"Total sequences in first file: {total1}")
    print(f"Total sequences in second file: {total2}")
    print(f"Number of common peptides: {len(common)}")
    print(f"Number of peptides unique to second FASTA: {len(only_second)}")

    # Write results to FASTA files
    with open("common_peptides.fasta", "w", encoding="utf-8") as out_c:
        for seq in sorted(common):
            header = seqs1.get(seq) or seqs2.get(seq)
            out_c.write(f"{header}\n{seq}\n")

    with open("unique_second.fasta", "w", encoding="utf-8") as out_u:
        for seq in sorted(only_second):
            header = seqs2[seq]
            out_u.write(f"{header}\n{seq}\n")

    print("Output files saved: common_peptides.fasta, unique_second.fasta")


if __name__ == "__main__":
    main()
