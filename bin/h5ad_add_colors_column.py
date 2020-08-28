#!/usr/bin/env python3

"""
Collaborators requested support for custom colors for their tSNE clusters but
we shouldn't have to re-generate datasets to supporr this.  This script takes
an input file of cluster names and appends a 'colors' column to the H5AD file.

"""

import argparse
import os
import pandas as pd
import scanpy as sc

def main():
    parser = argparse.ArgumentParser( description='Rename a colum in an H5AD file')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-ck', '--cluster_key', type=str, required=True, help='Name of column which identifies clusters')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_file)

    # https://stackoverflow.com/questions/26886653/pandas-create-new-column-based-on-values-from-other-columns-apply-a-function-o


    adata.write(args.output_file)

def choose_color(row):



if __name__ == '__main__':
    main()







