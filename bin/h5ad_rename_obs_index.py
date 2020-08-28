#!/usr/bin/env python3

"""
We commonly get a file uploaded to gEAR and converted to H5AD but there
was a mistake in its preparation and we want to rename an obs column without
redoing the entire process.

This script takes an H5AD and allows you to rename one column, writing a
new file.

Input: If your obs dataframe look slike this:

                        cell_type	louvain
index		
Outer_Pillar--1_	Outer_Pillar	Outer_Pillar
Outer_Pillar--2_	Outer_Pillar	Outer_Pillar
Outer_Pillar--3_	Outer_Pillar	Outer_Pillar
Outer_Pillar--4_	Outer_Pillar	Outer_Pillar
Outer_Pillar--5_	Outer_Pillar	Outer_Pillar

This script is for if you wanted to rename 'Outer_Pillar--5_' to 'Outer_Pillar_5', for 
example.
"""

import argparse
import os
import pandas as pd
import scanpy as sc

def main():
    parser = argparse.ArgumentParser( description='Rename a colum in an H5AD file')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-pre', '--previous_column', type=str, required=True, help='Name of column which needs changing')
    parser.add_argument('-post', '--post_column', type=str, required=True, help='New name for column')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_file)
    adata.obs.rename(index = {args.previous_column:args.post_column}, inplace = True)
    adata.write(args.output_file)


if __name__ == '__main__':
    main()







