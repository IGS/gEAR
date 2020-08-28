#!/usr/bin/env python3

"""

Returns the total number of columns and any rows that don't match that count

"""

import argparse
import csv
import json
import os
import re
import sys


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to file containing single cell data file' )
    args = parser.parse_args()
    filename = args.input_file

    f = open(filename)
    if filename.endswith('.tab'):
        print("Processing as tab-delimited file because of extension")
        reader = csv.reader(f, delimiter='\t')
    else:
        print("Processing as comma-delimited file because of extension")
        reader = csv.reader(f)

    line_num = 0
    col_count = 0

    for cols in reader:
        current_col_count = 0

        #Get Total number of columns
        if line_num == 0:
            col_count = len(cols)
            print("Total # of columns: {0}".format(col_count))
            line_num += 1
        else:

            # Get col count of current line
            current_col_count = len(cols)

            if current_col_count != col_count:
                gene = cols[0]
                print("{0} has {1} columns, but should have {2}".format(gene, current_col_count, col_count))

            line_num += 1


    f.close()

    print('\nDONE!')


if __name__ == '__main__':
    main()
