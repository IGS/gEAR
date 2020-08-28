#!/usr/bin/env python3

"""

This script is to convert an H5AD file into 3-tab output.

These three files are created in the output directory:

  expression.tab
  genes.tab
  observations.tab

INPUT

An HDF5 file in AnnData format

To convert to 10x:
https://github.com/theislab/scanpy/issues/262#issuecomment-549943099

pd.DataFrame(ad.var.index).to_csv(os.path.join(destination, "genes.tsv" ),   sep = "\t", index_col = False)
pd.DataFrame(ad.obs.index).to_csv(os.path.join(destination, "barcodes.tsv"), sep = "\t", index_col = False)
ad.obs.to_csv(os.path.join(destination, "metadata.tsv"), sep = "\t", index_col = True)
scipy.io.mmwrite(os.path.join(destination, "matrix.mtx"), ad.X)

"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from itertools import zip_longest
import os
import sys

def main():
    parser = argparse.ArgumentParser(description='H5AD to 3-tab output conversion script')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input h5ad file' )
    parser.add_argument('-o', '--output_directory', type=str, required=True, help='Path to an output directory in which to write the 3-tab files' )
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_file)
    adata.write_csvs(args.output_directory, sep="\t", skip_data=False)

    # https://github.com/theislab/scanpy/issues/506
    t=adata.X.toarray()
    pd.DataFrame(data=t, index=adata.obs_names, columns=adata.var_names).to_csv('adata_x.csv')

    sys.exit(0)

    print("adata.X is a {0}".format(type(adata.X)))

    # assumes adata.X is a <class 'scipy.sparse.csr.csr_matrix'>
    M = adata.X.tolil()

    matrix_tab = open("{0}/expression.tab".format(args.output_directory), 'wt')

    for cols in zip_longest(*M.data, fillvalue=''):
        matrix_tab.write("\t".join(str(v) for v in cols))

    sys.exit(0)

    # This parses/saves an HDF5 using scanpy
    adata = sc.read(args.input_directory + '/expression.tab', cache=False).transpose()
    obs = pd.read_table(args.input_directory + '/observations.tab', sep='\t', index_col=0, header=0)
    var = pd.read_table(args.input_directory + '/genes.tab', sep='\t', index_col=0, header=0)

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







