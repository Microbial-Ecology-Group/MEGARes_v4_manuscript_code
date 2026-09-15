#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def concatenate_csv_files(
    root_dir: str, output_file: str, include_source: bool = True
):
    """Recursively finds and concatenates all CSV files under root_dir.

    Args:
        root_dir (str): The parent directory containing subfolders with CSVs.
        output_file (str): The file path where the merged CSV will be saved.
        include_source (bool): If True, adds a column tracking the original file
        path.
    """
    root_path = Path(root_dir)
    output_path = Path(output_file)

    # 1. Recursively find all CSV files using a glob pattern
    csv_files = list(root_path.rglob("*.csv"))

    # Exclude the output file if it happens to be in the same path
    csv_files = [f for f in csv_files if f.resolve() != output_path.resolve()]

    if not csv_files:
        print(f"No CSV files found in '{root_dir}' or its subdirectories.")
        return

    print(f"📂 Found {len(csv_files)} CSV files. Commencing merge...")

    dfs = []
    for file_path in csv_files:
        try:
            # low_memory=False prevents mixed type warnings on large datasets
            df = pd.read_csv(file_path, low_memory=False)

            # Skip empty files
            if df.empty:
                continue

            # Optional: Add a metadata column to remember where the rows came from
            if include_source:
                df["source_file_path"] = str(
                    file_path.relative_to(root_path.parent)
                )

            dfs.append(df)
            print(f"✔️ Loaded: {file_path.name}")

        except Exception as e:
            print(f"⚠️ Failed to read {file_path.name}: {e}")

    # 2. Vertically stack all DataFrames
    if dfs:
        print("🔄 Concatenating data... Please wait.")
        # ignore_index=True re-indexes rows continuously from 0 to N
        combined_df = pd.concat(dfs, ignore_index=True, sort=False)

        # 3. Export the master CSV
        output_path.parent.mkdir(parents=True, exist_ok=True)
        combined_df.to_csv(output_path, index=False, encoding="utf-8-sig")
        print(f"Success! Merged data saved to: {output_path.resolve()}")
    else:
        print("No valid data could be retrieved from the discovered files.")


if __name__ == "__main__":
    # Setup Command Line Interface for maximum reusability
    parser = argparse.ArgumentParser(
        description="Concatenate CSV files from multiple subdirectories."
    )
    parser.add_argument(
        "--root",
        type=str,
        default=".",
        help="Path to the parent directory (default: current directory)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="merged_output.csv",
        help="Path/Name of the final output file",
    )
    parser.add_argument(
        "--keep-source",
        action="store_true",
        help="Include a column showing which file the row came from",
    )

    args = parser.parse_args()

    concatenate_csv_files(
        root_dir=args.root,
        output_file=args.output,
        include_source=args.keep_source,
    )

