#!/usr/bin/env python3

"""

This script is to convert three-tab text datafiles into an HDF5 AnnData file

These three files are expected to be found in the input directory:

  expression.tab
  genes.tab
  observations.tab

OUTPUT

An HDF5 file in AnnData format

If you only pass the --input_directory it assumes the files are called 'expression.tab',
'observations.tab' and 'genes.tab'

{input_directory}/expression.tab
{input_directory}/observations.tab
{input_directory}/genes.tab

If you also pass --file_prefix we're assuming NeMO file naming conventions and it will
instead look for:

{input_directory}/{file_prefix}_DataMTX.tab
{input_directory}/{file_prefix}_COLmeta.tab
{input_directory}/{file_prefix}_ROWmeta.tab
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import os


def main():
    parser = argparse.ArgumentParser(description='Creates an H5AD file from 3-tab input files')

    parser.add_argument('-i', '--input_directory', type=str, required=True, help='Path to an input directory containing the three tab files' )
    parser.add_argument('-f', '--file_prefix', type=str, required=False, help='Prefix portion for each of the three input files' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output H5AD file to be created' )
    args = parser.parse_args()

    if args.file_prefix:
        adata_file = "{0}/{1}_DataMTX.tab".format(args.input_directory, args.file_prefix)
        obs_file = "{0}/{1}_COLmeta.tab".format(args.input_directory, args.file_prefix)
        var_file = "{0}/{1}_ROWmeta.tab".format(args.input_directory, args.file_prefix)
    else:
        adata_file = args.input_directory + '/expression.tab'
        obs_file = args.input_directory + '/observations.tab'
        var_file = args.input_directory + '/genes.tab'

    # This parses/saves an HDF5 using scanpy
    adata = sc.read(adata_file, cache=False).transpose()
    obs = pd.read_table(obs_file, sep='\t', index_col=0, header=0)
    var = pd.read_table(var_file, sep='\t', index_col=0, header=0)

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







