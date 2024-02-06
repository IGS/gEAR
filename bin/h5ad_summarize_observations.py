#!/opt/bin/python3

"""
If you have an H5AD file with observation headers like this:

                                        cell_type  age                louvain
    observations                                                               
    AAACGCTAGCGAAACC_CA08_S1  Tympanic border cells  P12  Tympanic border cells
    AAGCGTTAGGAATCGC_CA08_S1             Osteocytes  P12             Osteocytes
    AATTTCCAGGGCCTCT_CA08_S1                B cells  P12                B cells
    ACGTAACAGACTAGAT_CA08_S1     Intermediate stria  P12     Intermediate stria
    ACTTTCATCATTTGGG_CA08_S1  Tympanic border cells  P12  Tympanic border cells
    AGCGTATAGGTAGTCG_CA08_S1  Tympanic border cells  P12  Tympanic border cells
    AGGGCCTGTCGCAACC_CA08_S1             Fibrocytes  P12             Fibrocytes
    AGGTTGTGTGGCGCTT_CA08_S1                    SS1  P12                    SS1
    AGTGTTGTCAGACAAA_CA08_S1  Tympanic border cells  P12  Tympanic border cells
    CACAGATCAATCAAGA_CA08_S1            Neutrophils  P12            Neutrophils

This script will print a summary of the observation keys you pass and the count of 
values in the column aggregate.  So if I passed keys: 'cell_type' and 'age', I would
get a summary like this:

    cell_type               age     count
    B cells                 P12     1
    Fibrocytes              P12     1
    Intermediate stria      P12     1
    Neutrophils             P12     1
    Osteocytes              P12     1
    SS1                     P12     1
    Tympanic border cells   P12     4
"""

import argparse
import os
import sys

from collections import defaultdict

import pandas as pd
import scanpy as sc

def main():
    parser = argparse.ArgumentParser( description='Summarize the observations (and count of each) in an H5AD file.')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-oc', '--obs_columns', type=str, required=True, help='Columns to group by and count, comma-separated')
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_file)

    passed_obs_columns = args.obs_columns.split(',')

    # make sure none of the passed observation columns are not in the adata.obs
    missing_columns = [c for c in passed_obs_columns if c not in adata.obs.columns]
    if missing_columns:
        print(f"Error: The following observation columns are not in the H5AD file: {missing_columns}")
        exit(1)

    # iterate over the adata.obs and count the values for each column
    grouping = adata.obs.groupby(passed_obs_columns).size()

    for i, v in grouping.iteritems():
        if isinstance(i, tuple):
            cols = list(i)
            cols.append(v)
        else:
            cols = [i, v]

        cols = [str(c) for c in cols]

        print("\t".join(cols))


if __name__ == '__main__':
    main()







