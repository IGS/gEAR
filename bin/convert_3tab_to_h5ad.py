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

def main():
    parser = argparse.ArgumentParser( description='3-tab -> H5AD')
    parser.add_argument('-i', '--input_directory', type=str, required=True, help='Path to an input directory containing the 3-tab files')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to the output file to be created')
    
    args = parser.parse_args()

    for infile in os.listdir(args.input_directory):
        filepath = "{0}/{1}".format(args.input_directory, infile)
        
        # Read each file as pandas dataframes
        if infile == 'expression.tab' or os.path.basename(filepath)== 'expression.tab' or 'DataMTX.tab' in infile:
            # Get columns and rows of expression data in list form.
            exp = pd.read_table(filepath, sep='\t', index_col=0, header=0)
            exp_obs = list(exp.columns)
            exp_genes= list(exp.index)
            exp = None

            # Read in expressions as AnnData object
            adata = sc.read(filepath, first_column_names=True, cache=False).transpose()
        elif infile == 'observations.tab' or os.path.basename(filepath)== 'observations.tab' or 'COLmeta.tab' in infile:
            obs = pd.read_table(filepath, sep='\t', index_col=0, header=0)
        elif infile == 'genes.tab' or os.path.basename(filepath)== 'genes.tab' or 'ROWmeta.tab' in infile:
            var = pd.read_table(filepath, sep='\t', index_col=0, header=0)

    for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
        if str_type in obs.columns:
            obs[str_type] = pd.Categorical(obs[str_type])

    for num_type in ['replicate', 'time_point_order']:
        if num_type in obs.columns:
            obs[num_type] = pd.to_numeric(obs[num_type])

    # Ensure observations and genes are sorted the same as found in the expressions file
    obs_index = list(obs.index)
    if set(obs_index) != set(exp_obs):
        print("Observation IDs from 'expressions' and 'observations' files are not the same.", file=sys.stderr)

        print("Got {0} obs IDs from observations file".format(set(obs_index)))
        print("Got {0} obs IDs from expressions file".format(set(exp_obs)))

        col_num = 0
    
        for (of_id, ef_id) in zip(set(obs_index), set(exp_obs)):
            col_num += 1

            if of_id != ef_id:
                print("Obs column mismatch at row/column {0}. Expression file has ({1}) but observation file has ({2})".format(col_num, ef_id, of_id))
                sys.exit(1)

        print("All observation values seemed to match")
        sys.exit(1)
    
    obs = obs.reindex(exp_obs)

    genes_index = list(var.index)
    if set(genes_index) != set(exp_genes):
        raise Exception("Gene IDs from 'expressions' and 'genes' files are not the same.")
    var = var.reindex(exp_genes)

    # Assign genes and observations to AnnData object
    adata.var = var
    adata.obs = obs

    adata.write(args.output_file)
            
if __name__ == '__main__':
    main()







