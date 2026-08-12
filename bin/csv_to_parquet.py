#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import argparse
import sys

"""
This script crowls through a specified directory and its subdirectories to find all CSV files,
converts them to Parquet format, and optionally deletes the original CSV files after successful conversion.
"""

def convert_directory(root_dir, delete_original=False):
    root_path = Path(root_dir)

    if not root_path.exists() or not root_path.is_dir():
        print(f"Error: Directory '{root_dir}' does not exist.")
        sys.exit(1)

    # Find all CSV files recursively
    csv_files = list(root_path.rglob("*.csv"))

    if not csv_files:
        print(f"No CSV files found in {root_dir}")
        return

    print(f"Found {len(csv_files)} CSV files. Starting conversion...\n")
    success_count = 0

    for csv_path in csv_files:
        try:
            # 1. Read the CSV
            df = pd.read_csv(csv_path)

            # 2. Datashader Optimization: Convert strings to Categorical
            # This automatically catches columns like 'cluster' or 'cell_type'
            for col in df.select_dtypes(include=['object']).columns:
                df[col] = df[col].astype('category')

            # 3. Define the new Parquet filepath
            parquet_path = csv_path.with_suffix('.parquet')

            # 4. Save as Parquet
            df.to_parquet(parquet_path, engine='pyarrow', index=False)

            # 5. Cleanup (Optional)
            if delete_original:
                csv_path.unlink()
                print(f"Converted & Deleted Original: {csv_path.name}")
            else:
                print(f"Converted: {csv_path.name} -> {parquet_path.name}")

            success_count += 1

        except Exception as e:
            print(f"Failed to convert {csv_path.name}: {e}")

    print(f"\nFinished! Successfully converted {success_count}/{len(csv_files)} files.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Recursively convert CSVs to Parquet files.")
    parser.add_argument("directory", help="Root directory to start crawling")
    parser.add_argument("--delete", action="store_true", help="Delete the original CSV files after successful conversion")

    args = parser.parse_args()

    # Run the converter
    convert_directory(args.directory, delete_original=args.delete)