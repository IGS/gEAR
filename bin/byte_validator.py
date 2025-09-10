#!/usr/bin/env python3

"""

A common problem is getting input files with non-UTF8 values which can kill the uploader.

We could modify to automatically replace them like this:

obs = pd.read_table(filepath, sep='\t', index_col=0, encoding='utf-8', errors='replace')

But that would insert a unicode replacement character like this 'badtext�more', so not 
sure if we want to do that.  This at least shows the rows, columns and actual bad values
so the user can replace them.

"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Reports position and values of non-UTF8 columns in an input file')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    args = parser.parse_args()

    bad_utf8 = validate_utf8_cells(args.input_file)
    

    

def validate_utf8_cells(filepath, max_report=100):
    """
    Scan a tab-delimited text file line by line and report cells that contain invalid UTF-8.
    Reports up to `max_report` bad cells.
    """
    bad_cells = []
    with open(filepath, 'rb') as f:
        for row_num, raw_line in enumerate(f, start=1):
            try:
                line = raw_line.decode('utf-8')
            except UnicodeDecodeError as e:
                # Whole line failed to decode — fallback to partial inspection
                line = None

            if line is not None:
                fields = line.rstrip('\n').split('\t')
                for col_num, field in enumerate(fields, start=1):
                    try:
                        field.encode('utf-8').decode('utf-8')  # Explicitly test
                    except UnicodeDecodeError:
                        bad_cells.append((row_num, col_num, repr(field)))
                        if len(bad_cells) >= max_report:
                            break
            else:
                # Decode failed — try recovering individual fields
                fields = raw_line.split(b'\t')
                for col_num, field in enumerate(fields, start=1):
                    try:
                        field.decode('utf-8')
                    except UnicodeDecodeError:
                        bad_cells.append((row_num, col_num, repr(field)))
                        if len(bad_cells) >= max_report:
                            break

            if len(bad_cells) >= max_report:
                break

    if bad_cells:
        print(f"Found {len(bad_cells)} bad UTF-8 cell(s):\n")
        for row, col, val in bad_cells:
            print(f"  🚨 Row {row}, Column {col}: {val}")
    else:
        print("✅ No bad UTF-8 found in file.")

    return bad_cells

if __name__ == '__main__':
    main()







