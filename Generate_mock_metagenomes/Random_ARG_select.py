#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 17:11:38 2026

@author: jared
"""

#!/usr/bin/env python3

import csv
import random
import argparse
from collections import defaultdict


def read_fasta(fasta_file):
    """
    Read a FASTA file into a dictionary:
    {header: sequence}

    Uses the first whitespace-delimited token in each FASTA header
    as the sequence ID, so:
      >gene123 some description
    becomes:
      gene123
    """
    sequences = {}
    header = None
    seq_chunks = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(seq_chunks)

                full_header = line[1:].strip()
                header = full_header.split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header is not None:
            sequences[header] = "".join(seq_chunks)

    return sequences


def read_annotations(annotation_file):
    """
    Read the annotation CSV and return rows as dictionaries.

    Expected columns:
      header, type, class, mechanism, group, snp
    """
    with open(annotation_file, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        expected = {"header", "type", "class", "mechanism", "group", "snp"}
        found = set(reader.fieldnames or [])

        if not expected.issubset(found):
            raise ValueError(
                f"CSV must contain columns {sorted(expected)}. "
                f"Found columns: {reader.fieldnames}"
            )

        rows = list(reader)

    return rows, reader.fieldnames


def filter_annotations(rows):
    """
    Exclude any rows where snp == 'RequiresSNPConfirmation'.
    Blank SNP values are allowed.
    """
    filtered = []
    for row in rows:
        snp_value = (row.get("snp") or "").strip()
        if snp_value == "RequiresSNPConfirmation":
            continue
        filtered.append(row)
    return filtered


def select_random_by_class(rows, n_per_class=10, seed=None):
    """
    Randomly select up to n_per_class rows for each unique class.
    """
    if seed is not None:
        random.seed(seed)

    grouped = defaultdict(list)
    for row in rows:
        class_name = (row.get("class") or "").strip()
        grouped[class_name].append(row)

    selected = []
    for class_name, class_rows in grouped.items():
        if len(class_rows) <= n_per_class:
            selected.extend(class_rows)
        else:
            selected.extend(random.sample(class_rows, n_per_class))

    return selected


def write_fasta(selected_rows, fasta_sequences, output_fasta):
    """
    Write selected sequences to a FASTA file.
    Returns:
      - found_headers: list of headers written
      - missing_headers: list of selected headers not found in FASTA
      - lengths: list of dicts with header and length
    """
    found_headers = []
    missing_headers = []
    lengths = []

    with open(output_fasta, "w") as out:
        for row in selected_rows:
            header = row["header"].strip()

            if header not in fasta_sequences:
                missing_headers.append(header)
                continue

            seq = fasta_sequences[header]
            found_headers.append(header)
            lengths.append({"header": header, "length": len(seq)})

            out.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

    return found_headers, missing_headers, lengths


def write_annotations(selected_rows, found_headers, fieldnames, output_csv):
    """
    Write only the selected annotation rows whose headers were found in the FASTA.
    """
    found_set = set(found_headers)

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for row in selected_rows:
            if row["header"].strip() in found_set:
                writer.writerow(row)


def write_lengths(lengths, output_csv):
    """
    Write header-length table.
    """
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["header", "length"])
        writer.writeheader()
        writer.writerows(lengths)


def main():
    parser = argparse.ArgumentParser(
        description="Randomly select up to 10 genes per class from annotations, "
                    "excluding RequiresSNPConfirmation, and extract matching FASTA entries."
    )
    parser.add_argument("fasta", help="Input database FASTA file")
    parser.add_argument("annotations", help="Input annotation CSV file")
    parser.add_argument(
        "--n-per-class",
        type=int,
        default=10,
        help="Number of genes to select per class (default: 10)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--out-fasta",
        default="limited_database.fasta",
        help="Output FASTA filename"
    )
    parser.add_argument(
        "--out-annotations",
        default="limited_annotations.csv",
        help="Output annotation CSV filename"
    )
    parser.add_argument(
        "--out-lengths",
        default="selected_gene_lengths.csv",
        help="Output gene lengths CSV filename"
    )

    args = parser.parse_args()

    fasta_sequences = read_fasta(args.fasta)
    annotation_rows, fieldnames = read_annotations(args.annotations)

    filtered_rows = filter_annotations(annotation_rows)
    selected_rows = select_random_by_class(
        filtered_rows,
        n_per_class=args.n_per_class,
        seed=args.seed
    )

    found_headers, missing_headers, lengths = write_fasta(
        selected_rows,
        fasta_sequences,
        args.out_fasta
    )

    write_annotations(
        selected_rows,
        found_headers,
        fieldnames,
        args.out_annotations
    )

    write_lengths(lengths, args.out_lengths)

    print(f"Total annotation rows read: {len(annotation_rows)}")
    print(f"Rows after SNP filter: {len(filtered_rows)}")
    print(f"Rows selected: {len(selected_rows)}")
    print(f"Sequences written to FASTA: {len(found_headers)}")
    print(f"Lengths written: {len(lengths)}")

    if missing_headers:
        print("\nWarning: The following selected headers were not found in the FASTA:")
        for h in missing_headers:
            print(f"  {h}")


if __name__ == "__main__":
    main()