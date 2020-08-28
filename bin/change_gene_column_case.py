#!/usr/bin/env python3

"""

Assumes a tab-delimited input file where the first column are gene symbols like this:

  Gene_symbol     Cancer--Cell_1  Cancer--Cell_10 Cancer--Cell_100
  AAGAB   154.561226827488        0.0     0.0
  AAR2    295.875190529996        299.455534712676        0.0
  AATF    546.792205537953        323.38381204192996      0.0

The --mode option can be used to tranform the case to 'upper', 'lower' or the default
mode of 'initcap', which would produce this:

  Gene_symbol     Cancer--Cell_1  Cancer--Cell_10 Cancer--Cell_100
  Aagab   154.561226827488        0.0     0.0
  Aar2    295.875190529996        299.455534712676        0.0
  Aatf    546.792205537953        323.38381204192996      0.0

"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## Add arguments for your script here.  Already added are -i and -o for input/output.
    #  More information on argument parsing here: https://docs.python.org/3/howto/argparse.html
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-m', '--mode', type=str, required=False, default='initcap', help='Mode of case in output: upper, lower, or initcap' )
    args = parser.parse_args()

    ofh = open(args.output_file, 'w')
    line_num = 0

    if args.mode not in ['upper', 'lower', 'initcap']:
        raise Exception("ERROR:  The supported modes are: 'upper', 'lower', 'initcap'")

    for line in open(args.input_file):
        line_num += 1
        
        line = line.rstrip()
        if len(line) == 0:
            continue

        cols = line.split("\t")

        if line_num > 1:
            if args.mode == 'initcap':
                cols[0] = cols[0].title()
            elif args.mode == 'lower':
                cols[0] = cols[0].lower()
            elif args.mode == 'upper':
                cols[0] = cols[0].upper()

        print("\t".join(cols), file=ofh)


if __name__ == '__main__':
    main()







