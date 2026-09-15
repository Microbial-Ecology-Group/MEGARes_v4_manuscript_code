#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 15:19:37 2026

@author: jared
"""

import argparse
import csv
import random
from pathlib import Path


def read_fasta_records(fasta_path):
    """
    Read a single-record or multi-record FASTA file.

    Returns a list of dictionaries:
        {
            "id": first word of FASTA header,
            "header": complete FASTA header without '>',
            "sequence": sequence
        }
    """
    records = []
    header = None
    sequence_parts = []

    with open(fasta_path, "r", encoding="utf-8-sig") as fasta_file:
        for line in fasta_file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append({
                        "id": header.split()[0],
                        "header": header,
                        "sequence": "".join(sequence_parts)
                    })

                header = line[1:].strip()
                sequence_parts = []

            else:
                sequence_parts.append(line.replace(" ", ""))

    if header is not None:
        records.append({
            "id": header.split()[0],
            "header": header,
            "sequence": "".join(sequence_parts)
        })

    if not records:
        raise ValueError(f"No FASTA records were found in {fasta_path}")

    ids = [record["id"] for record in records]

    if len(ids) != len(set(ids)):
        raise ValueError(
            f"Duplicate FASTA record IDs were found in {fasta_path}. "
            "The first word of each FASTA header must be unique."
        )

    return records


def write_fasta_records(records, output_path, line_width=60):
    """
    Write one or more FASTA records.
    """
    with open(output_path, "w") as output_file:
        for record in records:
            output_file.write(f">{record['header']}\n")

            sequence = record["sequence"]

            for index in range(0, len(sequence), line_width):
                output_file.write(
                    sequence[index:index + line_width] + "\n"
                )


def write_unused_genes(
    replacement_genes,
    assignments,
    output_directory
):
    """
    Write replacement genes that were not selected during the run.

    A gene is considered used if its unique FASTA record ID appears in at
    least one assignment. This also works with --allow_reuse: a gene selected
    multiple times is still counted as one used replacement gene.

    Returns:
        output_path: Path to the unused-gene FASTA
        unused_genes: list of unused FASTA records
    """
    used_gene_ids = {
        replacement_gene["id"]
        for site, replacement_gene in assignments
    }

    unused_genes = [
        gene
        for gene in replacement_genes
        if gene["id"] not in used_gene_ids
    ]

    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    output_filename = f"unused_genes_{len(unused_genes):03d}.fasta"
    output_path = output_directory / output_filename

    write_fasta_records(unused_genes, output_path)

    return output_path, unused_genes


def read_intervals_from_csv(
    csv_path,
    fragment_col,
    start_col,
    stop_col
):
    """
    Read removal intervals from a CSV.

    Coordinates are assumed to be 1-based and inclusive.

    Returns:
        {
            fragment_id: [(start, stop), ...]
        }
    """
    intervals_by_fragment = {}

    with open(
        csv_path,
        "r",
        newline="",
        encoding="utf-8-sig"
    ) as csv_file:

        reader = csv.DictReader(csv_file)

        if reader.fieldnames is None:
            raise ValueError("The coordinates CSV does not have a header row.")

        # Strip accidental whitespace from column names.
        reader.fieldnames = [
            column.strip() for column in reader.fieldnames
        ]

        required_columns = [
            fragment_col,
            start_col,
            stop_col
        ]

        missing_columns = [
            column
            for column in required_columns
            if column not in reader.fieldnames
        ]

        if missing_columns:
            raise ValueError(
                f"Missing CSV columns: {missing_columns}. "
                f"Found columns: {reader.fieldnames}"
            )

        for row_number, row in enumerate(reader, start=2):
            fragment_id = row[fragment_col].strip()

            if not fragment_id:
                raise ValueError(
                    f"Missing fragment ID in CSV row {row_number}."
                )

            try:
                start = int(row[start_col])
                stop = int(row[stop_col])
            except (ValueError, TypeError):
                raise ValueError(
                    f"Invalid start or stop coordinate in CSV row "
                    f"{row_number}: {row}"
                )

            if start < 1 or stop < 1:
                raise ValueError(
                    f"Coordinates must be at least 1 in CSV row "
                    f"{row_number}: {row}"
                )

            # Correct intervals entered in reverse order.
            if start > stop:
                start, stop = stop, start

            intervals_by_fragment.setdefault(
                fragment_id,
                []
            ).append((start, stop))

    return intervals_by_fragment


def merge_intervals(intervals):
    """
    Merge overlapping or directly adjacent intervals.

    Input and output coordinates are 1-based and inclusive.
    """
    if not intervals:
        return []

    sorted_intervals = sorted(
        intervals,
        key=lambda interval: interval[0]
    )

    merged = [sorted_intervals[0]]

    for current_start, current_stop in sorted_intervals[1:]:
        previous_start, previous_stop = merged[-1]

        if current_start <= previous_stop + 1:
            merged[-1] = (
                previous_start,
                max(previous_stop, current_stop)
            )
        else:
            merged.append((current_start, current_stop))

    return merged


def validate_intervals(intervals, sequence_length, fragment_id):
    """
    Confirm all intervals occur within the fragment.
    """
    for start, stop in intervals:
        if stop > sequence_length:
            raise ValueError(
                f"Interval ({start}, {stop}) exceeds the length of "
                f"fragment '{fragment_id}', which is "
                f"{sequence_length} bases."
            )


def assign_replacement_genes(
    insertion_sites,
    replacement_genes,
    random_generator,
    allow_reuse=False
):
    """
    Randomly assign replacement genes to insertion sites.

    Without reuse, every replacement gene can be selected only once.
    """
    number_of_sites = len(insertion_sites)
    number_of_genes = len(replacement_genes)

    if not allow_reuse and number_of_genes < number_of_sites:
        raise ValueError(
            f"There are {number_of_sites} insertion sites but only "
            f"{number_of_genes} replacement genes. Add more replacement "
            "genes or use --allow_reuse."
        )

    assignments = []

    if allow_reuse:
        for site in insertion_sites:
            gene = random_generator.choice(replacement_genes)
            assignments.append((site, gene))

    else:
        shuffled_genes = replacement_genes.copy()
        random_generator.shuffle(shuffled_genes)

        for site, gene in zip(
            insertion_sites,
            shuffled_genes[:number_of_sites]
        ):
            assignments.append((site, gene))

    return assignments


def replace_intervals_in_sequence(
    sequence,
    intervals_with_genes,
    genome_name,
    fragment_id
):
    """
    Replace intervals in one fragment with assigned genes.

    intervals_with_genes is a list containing:
        ((original_start, original_stop), replacement_gene)

    Returns:
        modified_sequence
        insertion_records
    """
    output_parts = []
    insertion_records = []

    # Position in the original sequence, using a 0-based Python index.
    original_cursor = 0

    # Number of bases currently present in the modified sequence.
    modified_length = 0

    for interval, replacement_gene in intervals_with_genes:
        original_start, original_stop = interval

        # Sequence before the interval being replaced.
        retained_sequence = sequence[
            original_cursor:original_start - 1
        ]

        output_parts.append(retained_sequence)
        modified_length += len(retained_sequence)

        replacement_sequence = replacement_gene["sequence"]

        # Final coordinates are 1-based and inclusive.
        inserted_start = modified_length + 1
        inserted_stop = (
            inserted_start + len(replacement_sequence) - 1
        )

        output_parts.append(replacement_sequence)
        modified_length += len(replacement_sequence)

        insertion_records.append({
            "genome": genome_name,
            "fragment": fragment_id,
            "original_removed_start": original_start,
            "original_removed_stop": original_stop,
            "original_removed_length": (
                original_stop - original_start + 1
            ),
            "inserted_gene_id": replacement_gene["id"],
            "inserted_gene_header": replacement_gene["header"],
            "inserted_gene_length": len(replacement_sequence),
            "inserted_start": inserted_start,
            "inserted_stop": inserted_stop,
            "strand": "+"
        })

        # Skip the original sequence through original_stop.
        original_cursor = original_stop

    # Add the sequence after the final replacement interval.
    trailing_sequence = sequence[original_cursor:]
    output_parts.append(trailing_sequence)

    return "".join(output_parts), insertion_records


def write_insertion_csv(insertion_records, output_csv):
    """
    Write insertion coordinates and replacement-gene information.
    """
    fieldnames = [
        "genome",
        "fragment",
        "original_removed_start",
        "original_removed_stop",
        "original_removed_length",
        "inserted_gene_id",
        "inserted_gene_header",
        "inserted_gene_length",
        "inserted_start",
        "inserted_stop",
        "strand"
    ]

    with open(
        output_csv,
        "w",
        newline="",
        encoding="utf-8"
    ) as csv_file:

        writer = csv.DictWriter(
            csv_file,
            fieldnames=fieldnames
        )

        writer.writeheader()
        writer.writerows(insertion_records)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Remove gene regions from a multi-fragment genome FASTA "
            "and replace them with randomly selected genes from a "
            "separate FASTA file."
        )
    )

    parser.add_argument(
        "genome_fasta",
        help="Input whole-genome FASTA, possibly containing multiple fragments"
    )

    parser.add_argument(
        "coordinates_csv",
        help="CSV containing fragment, start, and stop coordinates"
    )

    parser.add_argument(
        "replacement_genes_fasta",
        help="FASTA containing genes available for random insertion"
    )

    parser.add_argument(
        "output_fasta",
        help="Output FASTA containing the modified genome"
    )

    parser.add_argument(
        "output_csv",
        help="Output CSV containing inserted-gene coordinates"
    )

    parser.add_argument(
        "--fragment_col",
        default="fragment",
        help="CSV column containing the fragment or contig ID"
    )

    parser.add_argument(
        "--start_col",
        default="start",
        help="CSV column containing 1-based start coordinates"
    )

    parser.add_argument(
        "--stop_col",
        default="stop",
        help="CSV column containing 1-based inclusive stop coordinates"
    )

    parser.add_argument(
        "--genome_name",
        default=None,
        help=(
            "Genome name written to the output CSV. "
            "The input FASTA filename is used by default."
        )
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help=(
            "Random seed. Using the same seed reproduces the same "
            "replacement-gene assignments."
        )
    )

    parser.add_argument(
        "--allow_reuse",
        action="store_true",
        help=(
            "Allow the same replacement gene to be inserted at more "
            "than one site."
        )
    )

    parser.add_argument(
        "--header_suffix",
        default="_genes_replaced",
        help="Suffix added to each output FASTA header"
    )
    
    parser.add_argument(
        "--unused_genes_dir",
        default=None,
        help=(
            "Directory where the unused replacement-gene FASTA will be "
            "written. By default, it is written in the same directory as "
            "the output genome FASTA."
        )
    )

    args = parser.parse_args()

    genome_name = args.genome_name

    if genome_name is None:
        genome_name = Path(args.genome_fasta).stem

    random_generator = random.Random(args.seed)

    genome_records = read_fasta_records(args.genome_fasta)
    replacement_genes = read_fasta_records(
        args.replacement_genes_fasta
    )

    intervals_by_fragment = read_intervals_from_csv(
        args.coordinates_csv,
        args.fragment_col,
        args.start_col,
        args.stop_col
    )

    genome_ids = {
        record["id"] for record in genome_records
    }

    coordinate_fragment_ids = set(
        intervals_by_fragment.keys()
    )

    missing_fragments = coordinate_fragment_ids - genome_ids

    if missing_fragments:
        raise ValueError(
            "The following fragment IDs occur in the coordinates CSV "
            f"but not in the genome FASTA: {sorted(missing_fragments)}"
        )

    merged_intervals_by_fragment = {}
    insertion_sites = []

    for record in genome_records:
        fragment_id = record["id"]

        merged_intervals = merge_intervals(
            intervals_by_fragment.get(fragment_id, [])
        )

        validate_intervals(
            merged_intervals,
            len(record["sequence"]),
            fragment_id
        )

        merged_intervals_by_fragment[
            fragment_id
        ] = merged_intervals

        for interval in merged_intervals:
            insertion_sites.append({
                "fragment": fragment_id,
                "interval": interval
            })

    assignments = assign_replacement_genes(
        insertion_sites,
        replacement_genes,
        random_generator,
        allow_reuse=args.allow_reuse
    )
    
    if args.unused_genes_dir is None:
        unused_genes_directory = Path(args.output_fasta).parent
    else:
        unused_genes_directory = Path(args.unused_genes_dir)

    unused_genes_path, unused_genes = write_unused_genes(
        replacement_genes,
        assignments,
        unused_genes_directory
    )

    assignments_by_fragment = {}

    for site, replacement_gene in assignments:
        fragment_id = site["fragment"]
        interval = site["interval"]

        assignments_by_fragment.setdefault(
            fragment_id,
            []
        ).append((interval, replacement_gene))

    modified_records = []
    all_insertion_records = []

    total_original_length = 0
    total_modified_length = 0
    total_removed_length = 0
    total_inserted_length = 0

    for record in genome_records:
        fragment_id = record["id"]
        original_sequence = record["sequence"]

        fragment_assignments = assignments_by_fragment.get(
            fragment_id,
            []
        )

        # Ensure insertion sites are processed from left to right.
        fragment_assignments.sort(
            key=lambda assignment: assignment[0][0]
        )

        modified_sequence, insertion_records = (
            replace_intervals_in_sequence(
                original_sequence,
                fragment_assignments,
                genome_name,
                fragment_id
            )
        )

        modified_records.append({
            "id": fragment_id,
            "header": record["header"] + args.header_suffix,
            "sequence": modified_sequence
        })

        all_insertion_records.extend(insertion_records)

        total_original_length += len(original_sequence)
        total_modified_length += len(modified_sequence)

        for insertion_record in insertion_records:
            total_removed_length += insertion_record[
                "original_removed_length"
            ]
            total_inserted_length += insertion_record[
                "inserted_gene_length"
            ]

    write_fasta_records(
        modified_records,
        args.output_fasta
    )

    write_insertion_csv(
        all_insertion_records,
        args.output_csv
    )

    print(f"Genome name             : {genome_name}")
    print(f"Genome fragments        : {len(genome_records)}")
    print(f"Replacement genes       : {len(replacement_genes)}")
    print(f"Replacement sites       : {len(insertion_sites)}")
    print(f"Original genome length  : {total_original_length}")
    print(f"Removed sequence length : {total_removed_length}")
    print(f"Inserted sequence length: {total_inserted_length}")
    print(f"Modified genome length  : {total_modified_length}")
    print(f"Random seed             : {args.seed}")
    print(f"Output FASTA            : {args.output_fasta}")
    print(f"Output insertion CSV    : {args.output_csv}")
    print(f"Used replacement genes  : {len(replacement_genes) - len(unused_genes)}")
    print(f"Unused replacement genes: {len(unused_genes)}")
    print(f"Unused genes FASTA      : {unused_genes_path}")


if __name__ == "__main__":
    main()