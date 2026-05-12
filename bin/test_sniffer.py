#!/opt/bin/python3

"""
Test csv.Sniffer()

"""
import argparse
import csv
import os
import sys


def main():

    parser = argparse.ArgumentParser( description="Prints text STDOUT with csv.Sniffer()" )
    parser.add_argument('-f', '--input_file', type=str, required=False, help='Path to an input file (txt, tab, csv)' )
    args = parser.parse_args()


    if args.input_file is not None:
        print("ERROR: No input file included.\n")

    print("File: {0}".format(args.input_file))
    f = open(args.input_file)

    # detect delimiter before opening file
    # source: https://docs.python.org/3.5/library/csv.html#csv.Sniffer
    dialect = csv.Sniffer().sniff( f.read(1024) )
    f.seek(0)
    reader = csv.reader(f, delimiter=dialect.delimiter)

    line_count = 0
    for cols in reader:
        if line_count == 4:
            break
            
        print(cols)
        line_count += 1

    f.close()
    print("\nDone\n")

if __name__ == '__main__':
    main()
