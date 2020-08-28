#!/usr/bin/env python3

"""

Performs basic validation on an input TAB file.  Currently, this includes:

- Check that all rows have the same column count as the first row
- Check that all values outside the first row and column are integers or floats

"""

import argparse
import math
import os


def main():
    parser = argparse.ArgumentParser( description='Input tab file validator')

    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    args = parser.parse_args()

    first_col_count = None
    line_num = 1

    for line in open(args.input_file):
        cols = line.rstrip().split()

        if first_col_count is None:
            first_col_count = len(cols)
        else:
            if len(cols) != first_col_count:
                print("ERROR: {0} columns on line {1} instead of the expected {2}".format(len(cols), line_num, first_col_count))

        if line_num > 1:
            col_num = 1
            for col in cols[1:]:
                try:
                    new_num = float(col)
                    if math.isnan(new_num) or math.isinf(new_num):
                        print("ERROR: value in row {0}, column {1} doesn't appear to be numeric ({2})".format(line_num, col_num, col))
                except:
                    print("ERROR: value in row {0}, column {1} doesn't appear to be numeric ({2})".format(line_num, col_num, col))

                col_num += 1

        line_num += 1


if __name__ == '__main__':
    main()







