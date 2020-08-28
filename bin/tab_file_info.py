#!/usr/bin/env python3

"""
Quick utility to give stats on a passed tab file.  These include:

- Row count
- Min/Max column count
- Min/Max column length

"""

import sys

if len(sys.argv) != 2:
    print("\nUsage: {0} someinputfile.tab\n".format(sys.argv[0]), file=sys.stderr)
    exit(1)

row_count = 0
min_column_count=None
max_column_count=None
col_counts = dict()

min_column_length=None
max_column_length=None

for line in open(sys.argv[1]):
    row_count += 1

    cols = line.rstrip().split("\t")

    if min_column_count is None or min_column_count > len(cols):
        min_column_count = len(cols)

    if max_column_count is None or max_column_count < len(cols):
        max_column_count = len(cols)

    if len(cols) in col_counts:
        col_counts[len(cols)] += 1
    else:
        col_counts[len(cols)] = 1

    for cell in cols:
        if min_column_length is None or min_column_length > len(cell):
            min_column_length = len(cell)

        if max_column_length is None or max_column_length < len(cell):
            max_column_length = len(cell)


print("File: {0}\nRow_count: {1}\nMin/max column count: {2}/{3}".format(sys.argv[0], row_count, min_column_count,
                                                                          max_column_count))

if min_column_count != max_column_count:
    print("\tCols\tInstances:")
    for c in col_counts:
        print("\t{0}\t{1}".format(c, col_counts[c]))

print("Min/max column length: {0}/{1}".format(min_column_length,
                                              max_column_length
                                          ))
        
    
