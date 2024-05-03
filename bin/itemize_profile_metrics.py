#!/usr/bin/env python

"""
    itemize_profile_metrics.py

    This script computes information per dataset within a profile
"""

import sys
import argparse
from pathlib import Path

lib_path = Path(__file__).resolve().parents[1].joinpath('lib')
sys.path.append(str(lib_path))
import geardb


def main():
    parser = argparse.ArgumentParser(description="Compute the number of samples in a profile.")
    parser.add_argument("--share_id", "-s",  help="The share ID of the profile.")
    args = parser.parse_args()

    # Get the profile
    layout = geardb.get_layout_by_share_id(args.share_id)
    layout.get_members()

    profile_name = layout.label
    dataset_ids = layout.dataset_ids()

    print("num datasets:", len(dataset_ids))

    with open(f"{args.share_id}.tsv", "w") as f:
        f.write("Profile\tDataset\tSpecies\tDatatype\tNumGenes\tNumCells\n")
        dataset_collection = geardb.DatasetCollection()
        dataset_collection.get_by_dataset_ids(dataset_ids)

        for dataset in dataset_collection.datasets:
            print(f"Processing {dataset.id}")
            title = dataset.title
            title = f'"{title}"'
            species = dataset.organism
            datatype = dataset.dtype
            (num_samples, num_features) = dataset.get_shape(None, True)

            f.write(f"{profile_name}\t{title}\t{species}\t{datatype}\t{num_features}\t{num_samples}\n")

if __name__ == '__main__':
    main()