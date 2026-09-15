#!/usr/bin/env bash

set -euo pipefail

current_gene_pool="replacement_genes.fasta"
output_dir="replacement_results"

mkdir -p "$output_dir"

for genome_fasta in *.fasta; do
    genome_name="${genome_fasta%.fasta}"
    coordinates_csv="${genome_name}_start_stop_pos.csv"

    if [[ ! -f "$coordinates_csv" ]]; then
        echo "Skipping $genome_fasta: missing $coordinates_csv"
        continue
    fi

    genome_output_dir="${output_dir}/${genome_name}"
    mkdir -p "$genome_output_dir"

    echo "Processing $genome_name"
    echo "Current replacement pool: $current_gene_pool"

    python3 Whole_genome_cleaner_v2.py \
        "$genome_fasta" \
        "$coordinates_csv" \
        "$current_gene_pool" \
        "${genome_output_dir}/${genome_name}_replaced.fasta" \
        "${genome_output_dir}/${genome_name}_insertions.csv" \
        --genome_name "$genome_name" \
        --seed 42 \
        --unused_genes_dir "$genome_output_dir"

    unused_files=("${genome_output_dir}"/unused_genes_*.fasta)

    if [[ ${#unused_files[@]} -ne 1 || ! -f "${unused_files[0]}" ]]; then
        echo "Error: could not uniquely identify unused-gene FASTA for $genome_name"
        exit 1
    fi

    current_gene_pool="${unused_files[0]}"

    if [[ ! -s "$current_gene_pool" ]]; then
        echo "No replacement genes remain after $genome_name."
        break
    fi
done

echo "Final unused-gene pool: $current_gene_pool"
