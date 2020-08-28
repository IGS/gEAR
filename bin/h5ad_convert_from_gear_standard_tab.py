#!/usr/bin/env python3

"""

This script is to convert our common grouped format of gEAR data files into HDF5.  Headers on
the first line of the input would look like this:

   Ensembl_ID  Gene_symbol     Cancer--Cell_1  Cancer--Cell_10 Cancer--Cell_100   Mast--Cell_99   Mast--Cell_avg  myocyte--Cell_1 myocyte--Cell_10

Many columns are obviously omitted there, but the point is that it's in SITE--CELL_ID format.  If so, these
are parsed out.  If not,

OUTPUT

An HDF5 file in AnnData format

Also writes files under /tmp/


"""
import anndata
import argparse
import pandas as pd
import scanpy as sc
import os
import sys

def main():
    parser = argparse.ArgumentParser(description='Create H5AD from a tab file')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input TAB file' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output H5AD file to be created' )
    parser.add_argument('-d', '--dashed_columns', type=int, default=0, help='Set to 1 if columns are in SITE--CELL_ID format' )
    args = parser.parse_args()

    matrix_out = open("/tmp/expression.tab", 'wt')
    obs_out    = open("/tmp/observations.tab", 'wt')
    genes_out  = open("/tmp/genes.tab", 'wt')
    genes_out.write("Ensembl_ID\tgene_symbol\n")

    line_num = 0

    for line in open(args.input_file):
        line_num += 1
        line = line.rstrip()
        cols = line.split("\t")

        if len(cols) < 2: continue

        if line_num == 1:
            matrix_out.write("Ensembl_ID\t{0}\n".format("\t".join(cols[2:])))

            if args.dashed_columns:
                obs_out.write("observations\tcell_type\n")
                for col in cols[2:]:
                    try:
                        cond, label = col.split('--')
                    except ValueError:
                        print("Error on column:{0}".format(col))
                        sys.exit(1)

                    obs_out.write("{0}\t{1}\n".format(col, cond))
            else:
                obs_out.write("observations\n{0}".format("\n".join(cols[2:])))

        else:
            matrix_out.write("{0}\t{1}\n".format(cols[0], "\t".join(cols[2:])))
            try:
                genes_out.write("{0}\t{1}\n".format(cols[0], cols[1]))
            except IndexError:
                print("Error: this line didn't appear to have enough columns:\n{0}".format(line))

    matrix_out.close()
    obs_out.close()
    genes_out.close()

    foo = input("Modify any of the files in /tmp as needed, then hit <enter> to finish")

    exp_df = pd.read_table("/tmp/expression.tab", index_col=0).transpose()
    obs_df = pd.read_table("/tmp/observations.tab", index_col=0)
    genes_df = pd.read_table("/tmp/genes.tab", index_col=0)

    X = exp_df.values[:, 0:].astype(float)
    adata = anndata.AnnData(X=X, obs=obs_df, var=genes_df)
    adata.write(filename=args.output_file)


if __name__ == '__main__':
    main()







