#!/usr/bin/env python3

"""
We commonly get a file uploaded to gEAR and converted to H5AD but there
was a mistake in its preparation and we want to perform a find/replace
in the values of a column.

This script takes an H5AD and allows you to rename all values in an obs
column.  The -i and -o can be the same file.


"""

import argparse
import os
import pandas as pd
import scanpy as sc

def main():
    parser = argparse.ArgumentParser( description='Rename a colum in an H5AD file')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-c', '--column', type=str, required=True, help='Column name in obs to find/search' )
    parser.add_argument('-pre', '--previous_value', type=str, required=True, help='Current value to search for')
    parser.add_argument('-post', '--post_value', type=str, required=True, help='New value to replace with')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_file)

    if type(adata.obs[args.column]):
        adata.obs[args.column].replace({args.previous_value:args.post_value}, inplace=True)
        adata.write(args.output_file)
    else:
        print("ERROR: Didn't find column '{0}' in the dataset observations".format(args.column))

if __name__ == '__main__':
    main()







