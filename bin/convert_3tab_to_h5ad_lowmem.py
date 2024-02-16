#!/opt/bin/python3

"""
This script takes as input a directory which contains these three files:

  expression.tab
  genes.tab
  observations.tab

And does the same process as the web-based uploader to convert it to an H5AD file.

"""

import argparse
import os
import re
import sys

import pandas as pd
import scanpy as sc
from scipy import sparse


def main():
    parser = argparse.ArgumentParser( description='3-tab -> H5AD')
    parser.add_argument('-i', '--input_directory', type=str, required=True, help='Path to an input directory containing the 3-tab files')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to the output file to be created')
    
    args = parser.parse_args()

    expression_matrix_path = None
    obs = None
    var = None

    for infile in os.listdir(args.input_directory):
        filepath = "{0}/{1}".format(args.input_directory, infile)
        
        # Read each file as pandas dataframes
        if infile == 'expression.tab' or os.path.basename(filepath)== 'expression.tab' or 'DataMTX.tab' in infile:
            expression_matrix_path = filepath
        elif infile == 'observations.tab' or os.path.basename(filepath)== 'observations.tab' or 'COLmeta.tab' in infile:
            print("Reading observations file: {0}".format(filepath), file=sys.stderr, flush=True)
            obs = pd.read_table(filepath, sep='\t', index_col=0, header=0)
        elif infile == 'genes.tab' or os.path.basename(filepath)== 'genes.tab' or 'ROWmeta.tab' in infile:
            print("Reading genes file: {0}".format(filepath), file=sys.stderr, flush=True)
            var = pd.read_table(filepath, sep='\t', index_col=0, header=0)

    # Read in expressions as AnnData object in a memory-efficient manner
    print("Creating AnnData object with obs and var", file=sys.stderr, flush=True)
    adata = sc.AnnData(obs=var, var=obs)
    print("Reading expression matrix file: {0}".format(expression_matrix_path), file=sys.stderr, flush=True)    
    reader = pd.read_csv(expression_matrix_path, sep='\t', index_col=0, chunksize=500)
    adata.X = sparse.vstack([sparse.csr_matrix(chunk.values) for chunk in reader])
    print("Finished reading expression matrix file", file=sys.stderr, flush=True)

    for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
        if str_type in obs.columns:
            obs[str_type] = pd.Categorical(obs[str_type])

    for num_type in ['replicate', 'time_point_order']:
        if num_type in obs.columns:
            obs[num_type] = pd.to_numeric(obs[num_type])

    # Assign genes and observations to AnnData object
    adata.var = var
    adata.obs = obs

    adata.write(args.output_file)
            
if __name__ == '__main__':
    main()







