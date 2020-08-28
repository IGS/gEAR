#!/opt/bin/python3

"""
This takes an input tab file and creates a three-file MEX for import into gEAR.

Input tab file should have columns:

1. Ensembl_id
2. Gene symbol
3-...  Expression columns

"""

import argparse
import os

def main():
    parser = argparse.ArgumentParser( description='Create MEX files from gEAR standard tab')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_directory', type=str, required=True, help='Path to an output directory where MEX files will be written' )
    args = parser.parse_args()

    barcodes_fh = open("{0}/barcodes.tsv".format(args.output_directory), 'w')
    genes_fh = open("{0}/genes.tsv".format(args.output_directory), 'w')
    matrix_fh = open("{0}/matrix.mtx".format(args.output_directory), 'w')

    line_number = 0

    for line in open(args.input_file):
        line_number += 1
        line = line.rstrip()
        cols = line.split("\t")

        if line_number == 1:
            barcodes_fh.write("{0}\n".format("\n".join(cols[2:])))

        genes_fh.write("{0}\t{1}\n".format(cols[0], cols[1]))
        matrix_fh.write("{0}\t{1}\n".format(cols[0], "\t".join(cols[2:])))

if __name__ == '__main__':
    main()







