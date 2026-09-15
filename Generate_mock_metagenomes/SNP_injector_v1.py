#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:00:50 2026

@author: jared
"""

#!/usr/bin/env python3

import argparse
import csv
import random
from typing import List, Tuple


VALID_BASES = {"A", "C", "G", "T"}


def read_fasta(fasta_path: str) -> List[Tuple[str, str]]:
    """
    Read a FASTA file and return a list of (header, sequence) tuples.
    The header is stored without the leading '>'.
    """
    records = []
    header = None
    seq_chunks = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header is not None:
            records.append((header, "".join(seq_chunks)))

    return records


def write_fasta(records: List[Tuple[str, str]], output_path: str, line_width: int = 80) -> None:
    """
    Write FASTA records to file.
    """
    with open(output_path, "w") as f:
        for header, sequence in records:
            f.write(f">{header}\n")
            for i in range(0, len(sequence), line_width):
                f.write(sequence[i:i + line_width] + "\n")


def choose_alt_base(ref_base: str) -> str:
    """
    Choose a nucleotide different from ref_base.
    """
    ref_base = ref_base.upper()
    choices = [b for b in ["A", "C", "G", "T"] if b != ref_base]
    return random.choice(choices)


def inject_snps_into_sequence(
    header: str,
    sequence: str,
    min_snps: int = 1,
    max_snps: int = 10
) -> Tuple[str, List[dict]]:
    """
    Inject random SNPs into a single sequence.

    Returns:
        mutated_sequence
        list of SNP records with 0-based positions
    """
    seq_list = list(sequence.upper())

    # Only allow SNP injection at standard nucleotide positions
    candidate_positions = [i for i, base in enumerate(seq_list) if base in VALID_BASES]

    if not candidate_positions:
        # No valid A/C/G/T positions to mutate
        return "".join(seq_list), []

    max_possible = min(max_snps, len(candidate_positions))
    min_possible = min(min_snps, max_possible)

    if max_possible == 0:
        return "".join(seq_list), []

    n_snps = random.randint(min_possible, max_possible)
    chosen_positions = random.sample(candidate_positions, n_snps)

    snp_records = []

    for pos in sorted(chosen_positions):
        ref_base = seq_list[pos]
        alt_base = choose_alt_base(ref_base)
        seq_list[pos] = alt_base

        snp_records.append({
            "header": header,
            "reference_nucleotide": ref_base,
            "zero_based_position": pos,
            "one_based_position": pos + 1,
            "variant_nucleotide": alt_base,
        })

    return "".join(seq_list), snp_records


def write_zero_based_csv(records: List[dict], output_path: str) -> None:
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["header", "reference_nucleotide", "zero_based_position", "variant_nucleotide"])
        for r in records:
            writer.writerow([
                r["header"],
                r["reference_nucleotide"],
                r["zero_based_position"],
                r["variant_nucleotide"]
            ])


def write_one_based_csv(records: List[dict], output_path: str) -> None:
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["header", "reference_nucleotide", "one_based_position", "variant_nucleotide"])
        for r in records:
            writer.writerow([
                r["header"],
                r["reference_nucleotide"],
                r["one_based_position"],
                r["variant_nucleotide"]
            ])


def main():
    parser = argparse.ArgumentParser(
        description="Inject 1-10 random artificial SNPs into each gene in a FASTA file."
    )
    parser.add_argument(
        "input_fasta",
        help="Input FASTA file containing gene sequences"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for reproducibility"
    )
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    input_records = read_fasta(args.input_fasta)

    mutated_records = []
    all_snp_records = []

    for header, sequence in input_records:
        mutated_seq, snp_records = inject_snps_into_sequence(header, sequence)
        mutated_records.append((header, mutated_seq))
        all_snp_records.extend(snp_records)

    write_fasta(mutated_records, "SNP_injected_genes.fasta")
    write_zero_based_csv(all_snp_records, "zero_based_SNP_key.csv")
    write_one_based_csv(all_snp_records, "one_based_SNP_key.csv")

    print(f"Processed {len(input_records)} gene(s).")
    print("Created:")
    print("  SNP_injected_genes.fasta")
    print("  zero_based_SNP_key.csv")
    print("  one_based_SNP_key.csv")


if __name__ == "__main__":
    main()