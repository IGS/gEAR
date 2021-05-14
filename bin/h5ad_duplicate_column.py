#!/usr/bin/env python3

"""
This script takes an H5AD and allows you to duplicate a column's values but give
it a new name.
"""

import argparse
import os
import pandas as pd
import scanpy as sc

def main():
    parser = argparse.ArgumentParser( description='Rename a colum in an H5AD file')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-s', '--source', type=str, required=True, help='Source column name to copy')
    parser.add_argument('-d', '--destination', type=str, required=True, help='Destination column name')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_file)

    adata.obs[args.destination] = adata.obs[args.source]
    adata.write(args.output_file)

if __name__ == '__main__':
    main()







