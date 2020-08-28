#!/usr/bin/env python3

"""
Prints out a tab-delimited report of all datasets and their sizes.  Columns:

Dataset_ID      row_count       column_count      cell_count

The last line says "SUMMARY" for the ID and is the sum of all datasets.

Note, this does NOT check for empty columns, only counts the overall dimensions.
"""

import scanpy as sc
import pandas as pd
import os
import sys

dataset_dir = sys.argv[1]

if not dataset_dir:
    print("\nUsage: report_dataset_dimensions.py [directory with h5ad files]")
    sys.exit(1)

sum_rows = 0
sum_cols = 0
sum_cells = 0

for file in os.listdir(dataset_dir):
    if file.endswith('.h5ad'):
        filepath = "{0}/{1}".format(dataset_dir, file)
        dataset_id = file.rstrip('.h5ad')
        adata = sc.read_h5ad(filepath, backed='r')
        cols_X, rows_X = adata.X.shape
        cells_X = cols_X * rows_X
        sum_rows += rows_X
        sum_cols += cols_X
        sum_cells += cells_X
        print("{0}\t{1}\t{2}\t{3}".format(dataset_id, rows_X, cols_X, cells_X))

print("{0}\t{1}\t{2}\t{3}".format('SUMMARY', sum_rows, sum_cols, sum_cells))


