#!/opt/bin/python3

"""
If the input tab file has unconventional empty columns (such as the "NA" string)
this can be used to empty it or set it to something else.

Skips first/header line
"""

import argparse
import os, sys

lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'lib')
sys.path.append(lib_path)
import geardb


def main():
    parser = argparse.ArgumentParser( description='Resets column values matching one thing to another value')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-is', '--input_string', type=str, required=True, help='String to look for (full column values only)')
    parser.add_argument('-os', '--output_string', type=str, required=True, help='String to write instead')
    args = parser.parse_args()

    ofh = open(args.output_file, 'wt')
    
    line_count = 0
    
    for line in open(args.input_file):
        line_count += 1
        line = line.rstrip()
        
        if line_count == 1:
            print(line, file=ofh)
            continue

        cols = list()
        
        for col in line.split("\t"):
            if col == args.input_string:
                cols.append(args.output_string)
            else:
                cols.append(col)

        print("\t".join(cols), file=ofh)
        

if __name__ == '__main__':
    main()
