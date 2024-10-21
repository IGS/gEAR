#!/usr/bin/env python3

"""
Print best match ensembl release for set of ensembl ids from file, count, and organism id.

Assumes Ensembl ID or gene symbol is the first column of the input file
"""

import argparse
import os, sys

lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'lib')
sys.path.append(lib_path)
from gear.orthology import find_best_ensembl_release_match

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read.')
    parser.add_argument('-org', '--organism_id', type=int, required=True, help='Organism ID')
    parser.add_argument("--silent", action="store_true", required=False, help="Suppress printing. Only return best release version")
    args = parser.parse_args()
    find_best_ensembl_release_match(args.input_file, args.organism_id, args.silent)

if __name__ == '__main__':
    main()