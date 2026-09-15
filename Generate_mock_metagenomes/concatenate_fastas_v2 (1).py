#!/usr/bin/env python3

from pathlib import Path
import argparse


def concatenate_fasta_files(parent_directory, output_file):
    """
    Recursively find and concatenate FASTA files from all subdirectories
    within a parent directory.
    """

    parent_directory = Path(parent_directory)
    output_file = Path(output_file).resolve()

    if not parent_directory.is_dir():
        raise NotADirectoryError(
            f"Directory does not exist: {parent_directory}"
        )

    fasta_files = sorted(
    	file
    	for file in parent_directory.rglob("*_replaced.fasta")
    	if file.resolve() != output_file
    )

    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files were found beneath: {parent_directory}"
        )

    with output_file.open("w") as outfile:
        for fasta_file in fasta_files:
            print(f"Adding: {fasta_file}")

            with fasta_file.open("r") as infile:
                contents = infile.read()

                outfile.write(contents)

                # Ensure the next FASTA file begins on a new line.
                if contents and not contents.endswith("\n"):
                    outfile.write("\n")

    print(
        f"\nCombined {len(fasta_files)} FASTA files "
        f"into: {output_file}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Recursively concatenate FASTA files from subdirectories "
            "into a single FASTA file."
        )
    )

    parser.add_argument(
        "parent_directory",
        help="Parent directory containing the FASTA subdirectories"
    )

    parser.add_argument(
        "output_file",
        help="Path for the combined FASTA output file"
    )

    args = parser.parse_args()

    concatenate_fasta_files(
        args.parent_directory,
        args.output_file
    )
