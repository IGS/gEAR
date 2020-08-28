#!/usr/bin/env python3

"""

Search the dataset file for rows containing the word given, ex 'NA' or 'N/A',
and remove that row.

Meant for gEAR datasets prior to being uploaded to the gEAR Portal.

File MUST be in TAB FORMAT.

"""

import argparse
import os
import json
import sys
import re


def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    # output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input TAB file to be read' )
    parser.add_argument('-o', '--output_file_name', type=str, required=True, help='Name of output file' )
    parser.add_argument('-r', '--remove_word', type=str, required=True, help='Rows containing this word will be removed, ex: ' )
    args = parser.parse_args()


    pattern = "\t" + args.remove_word 
    filename = args.input_file
    output_name = "./{0}".format(args.output_file_name)
    

    encoding = 'utf-8'
    matched = re.compile(pattern).search
    rows_deleted = 0

    with open(filename, encoding=encoding) as input_file:
        output_file = open(output_name, "w")
        for line in input_file:
            if not matched(line):
                output_file.write(line)
            if matched(line):
                rows_deleted += 1
    print('\nDONE!\n{0} rows were found to contain "{1}" and were removed.'.format(rows_deleted, args.remove_word))


if __name__ == '__main__':
    main()







