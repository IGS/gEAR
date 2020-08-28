#!/usr/bin/env python3

"""

This script is to convert Market Exchange Format (MEX) datafiles into an HDF5 AnnData file

These three files are expected to be found in the input directory:

  matrix.mtx
  genes.tsv
  barcodes.tsv

OUTPUT

An HDF5 file in AnnData format


"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import os


def main():
    parser = argparse.ArgumentParser(description='Put a description of your script here')

    parser.add_argument('-i', '--input_directory', type=str, required=True, help='Path to an input TAB file' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output H5AD file to be created' )
    args = parser.parse_args()

    # This parses/saves an HDF5 using scanpy
    adata = sc.read(args.input_directory + '/matrix.mtx', cache=True).transpose()
    adata.var_names = np.genfromtxt(args.input_directory + '/genes.tsv', dtype=str)[:, 1]
    adata.obs_names = np.genfromtxt(args.input_directory + '/barcodes.tsv', dtype=str)

    adata.write(args.output_file)


if __name__ == '__main__':
    main()







