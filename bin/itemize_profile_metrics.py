#!/usr/bin/env python

"""
    itemize_profile_metrics.py

    This script computes information per dataset within a profile
"""

import sys
import argparse
from pathlib import Path

lib_path = Path(__file__).resolve().parents[1].joinpath('lib')
sys.path.append(lib_path)
import geardb


def main():
    parser = argparse.ArgumentParser(description="Compute the number of samples in a profile.")
    parser.add_argument("share_id", "-s",  help="The share ID of the profile.")
    args = parser.parse_args()

    # Get the profile
    layout = geardb.get_layout_by_share_id(args.share_id)
    layout.get_members()

    profile_name = layout.label

    if args.unique_datasets:
        dataset_ids = layout.dataset_ids()

    with open(f"{args.share_id}.tsv", "w") as f:
        f.write("Profile\tDataset\tSpecies\tDatatype\tNumGenes\tNumCells\t")
        dataset_collection = DatasetCollection()
        dataset_collection.get_datasets_by_ids(dataset_ids)

        for dataset in dataset_collection.datasets:
            species = dataset.organism
            dtype = dataset.datatype
            (num_samples, num_features) = dataset.get_shape(None, True)

            f.write(f"{profile_name}\t{dataset.name}\t{species}\t{dtype}\t{num_features}\t{num_samples}\n")

if __name__ == '__main__':
    main()