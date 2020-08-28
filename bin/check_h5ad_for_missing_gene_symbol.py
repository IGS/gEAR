#!/usr/bin/env python3

"""
Given a path to a directory with h5ad files, check if h5ads
have gene_symbol column in .var
"""

from glob import glob
import scanpy as sc
import argparse

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-d',
    '--directory',
    type=str,
    required=True,
    help='Path to directory of H5ads.' )
  args = parser.parse_args()

  h5ads = glob(f'{args.directory}/*.h5ad')
  for h5ad in h5ads:
    adata = sc.read(h5ad)
    if 'gene_symbol' in adata.var.columns:
      print(f'{h5ad} ✅')
    else:
      print(f'{h5ad} ❌')

if __name__ == "__main__":
    main()
