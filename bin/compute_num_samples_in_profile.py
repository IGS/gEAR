#!/usr/bin/env python

"""
    compute_num_samples_in_profile.py

    This script computes the number of samples of all datasets in a profile.
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
    parser.add_argument("unique_datasets", "-u", action="store_true", default=False, help="Only count unique datasets.")
    args = parser.parse_args()

    # Get the profile
    layout = geardb.get_layout_by_share_id(args.share_id)
    layout.get_members()
    if args.unique_datasets:
        dataset_ids = layout.dataset_ids()
    else:
        members = layout.members
        dataset_ids = [member.dataset_id for member in members]

    # Compute the number of samples
    total_samples = 0
    for dataset_id in dataset_ids:
        dataset = geardb.get_dataset_by_id(dataset_id)
        (num_samples, num_features) = dataset.get_shape(None, True)
        total_samples += num_samples

    print(layout.label)
    print(total_samples)

if __name__ == '__main__':
    main()