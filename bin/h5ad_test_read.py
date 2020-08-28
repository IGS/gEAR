#!/usr/bin/env python3

"""
Prints out general structural information about a passed H5AD file.
"""

import scanpy as sc
import pandas as pd
import os
import sys

h5_files_to_search = list()

filename = sys.argv[1]

if len(sys.argv) == 3:
    print_col = sys.argv[2]
else:
    print_col = None

if not filename:
    print("\nUsage: h5ad_test_read.py /path/to/some_file.h5ad [column_name_to_print]")
    sys.exit(1)

adata = sc.read_h5ad(filename, backed='r')

cols_X, rows_X = adata.X.shape
print("X:\n\tshape: {0} rows, {1} columns".format(rows_X, cols_X))

print("\n---------------------------------------------")
print("---------------------------------------------")

print("Observation columns:")
for col in adata.obs.columns:
    print("\t{0}".format(col))

print(adata.obs.head(5))

print("\n---------------------------------------------")
print("---------------------------------------------")

print("Var columns:")
for col in adata.var.columns:
    print("\t{0}".format(col))

# if you want to print all genes
#with pd.option_context('display.max_rows', None, 'display.max_columns', None):
#    print(adata.var)

print(adata.var.head(5))

print("\n---------------------------------------------")
print("------------------ Analyses -----------------")
print("---------------------------------------------\n")

analysis_count = 0

try:
    if type(adata.obsm['X_tsne']):
        print("tSNE")
        analysis_count += 1

except:
    pass

if not analysis_count:
    print("\tNone found")

if print_col:
    print("\n-------------------------------------------------")
    print("------------------ Values in '{0}' -----------------".format(print_col))
    print("-------------------------------------------------\n")

    if print_col in adata.obs:
        col_vals = set(adata.obs[print_col].tolist())

        for val in sorted(col_vals):
            print("\t{0}".format(val))
    else:
        print("\nNot found")
