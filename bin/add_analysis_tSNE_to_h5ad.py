#!/usr/bin/env python3

"""

Assumes the input tab file has columns named tSNE1 and tSNE2.  The first column of
the file should be the column index matching the data matrix.

H5AD storage:
adata.obsm['X_tsne'] =

array([[-13.04915115, -20.2436883 ],
       [  3.30309571,   8.50967076],
       [  2.73049336,   6.58619275],
       [ -0.19263491,  14.53663752],
       ...
       [  4.86175086,   5.51422638]])
"""

import argparse
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import os
import sys

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')
    parser.add_argument('-i', '--input_h5ad', type=str, required=True, help='Path to an input H5AD file to be read' )
    parser.add_argument('-t', '--input_tsne', type=str, required=True, help='Input tab-delimited file with columns: tSNE1, tSNE2' )
    parser.add_argument('-o', '--output_h5ad', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    # read the input file, find the tSNE1 and tSNE2 columns, get keyed arrays
    tsne_idx = get_tsne_idx(args.input_tsne)

    adata = sc.read_h5ad(args.input_h5ad)
    tsne_coords = list()

    for cell_name in adata.obs.index:
        if cell_name not in tsne_idx:
            raise Exception("Error: cell name ({}) in H5AD wasn't found in the tab input index column".format(cell_name))

        tsne_coords.append(tsne_idx[cell_name])

    adata.obsm['X_tsne'] = ad.base.BoundRecArr(input_array=np.array(tsne_coords), parent=adata, attr='obsm')

    adata.write(args.output_h5ad)

def get_tsne_idx(ifile):
    idx = dict()
    line_num = 0
    tSNE1_col = None
    tSNE2_col = None

    for line in open(ifile):
        line = line.rstrip()
        cols = line.split("\t")

        if line_num == 0:
            col_num = 0
            for col in cols:
                if col == 'tSNE1':
                    tSNE1_col = col_num
                elif col == 'tSNE2':
                    tSNE2_col = col_num

                col_num += 1

            # we should have found both
            if tSNE1_col is None:
                raise Exception("Error: failed to find 'tSNE1' column on first line of --input_tsne file")

            if tSNE2_col is None:
                raise Exception("Error: failed to find 'tSNE2' column on first line of --input_tsne file")
        else:
            idx[cols[0]] = [float(cols[tSNE1_col]), float(cols[tSNE2_col])]

        line_num += 1

    return idx

if __name__ == '__main__':
    main()
