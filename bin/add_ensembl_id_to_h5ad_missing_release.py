#!/opt/bin/python3

"""
Given an h5ad file with only gene symbols, find best matching ensembl release, then
collect ensembl ids for the corresponding gene symbol.

Fake ensembl IDs will be created for those genes whose gene symbols don't map to
a known ensembl ID.
"""

import argparse
import scanpy as sc
import sys
import os

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)

from utils import update_adata_with_ensembl_ids


def main():
    parser = argparse.ArgumentParser(
        description='See the script name?  It does that.')
    parser.add_argument('-i', '--input_file', type=str,
                        required=True, help='Path to an h5ad file to be read.')
    parser.add_argument('-o', '--output_file', type=str,
                        required=True, help='Output file to be written.')
    parser.add_argument('-org', '--organism', type=int,
                        required=True, help='Organism ID.')
    parser.add_argument('-r', '--read_only', type=bool, default=False,
                        help='If true, only print and do not write to output.')
    parser.add_argument('-idp', '--id_prefix', type=str,
                        required=True, help="Fake Ensembl IDs will be generated for any rows whose gene symbols don't match, using this prefix.")

    args = parser.parse_args()

    adata = sc.read(args.input_file)

    print('############################################################')
    print(f"File: {args.input_file}")
    adata = update_adata_with_ensembl_ids(adata, args.organism, args.id_prefix)

    print("FINAL ADATA")
    print(adata)
    print('VAR\n')
    print(adata.var)
    print("OBS\n")
    print(adata.obs)


    if not args.read_only:
        adata.write(args.output_file)
    print('############################################################')


if __name__ == '__main__':
    main()
