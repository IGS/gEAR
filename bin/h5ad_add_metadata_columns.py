#!/usr/bin/env python3

"""

"""

import argparse
import os
import pandas as pd
import scanpy as sc

def main():
    parser = argparse.ArgumentParser( description='Rename a colum in an H5AD file')
    parser.add_argument('-im', '--input_matrix', type=str, required=True, help='Path to an input h5ad file to be read' )
    parser.add_argument('-mf', '--metadata_file', type=str, required=True, help='Path to input file of new metadata columns to be added. Index should correspond to adata.obs.index')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input_matrix)

    # read a tab-delimited input file, set index to first column
    metadata = pd.read_csv(args.metadata_file, sep='\t')
    metadata = metadata.set_index(metadata.columns[0])

    # make sure indexes match. If they don't, show where they disagree
    assert (adata.obs.index == metadata.index).all()

    # See where the indexes disagree
    #print(adata.obs.index.difference(metadata.index))

    # add metadata columns to adata
    adata.obs = adata.obs.join(metadata)

    # remove column of a certain name
    #adata.obs.drop('Unnamed..0', axis=1, inplace=True)

    adata.write(args.output_file)


if __name__ == '__main__':
    main()







