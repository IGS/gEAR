#!/usr/bin/env python3

"""

Scanpy doesn't play well with input matrices which have numeric headers.

This script is used to read in these and replace them with string values by
prepending text to them.

These three files are expected to be found in the input directory:

  expression.tab
  genes.tab
  observations.tab

Or their NeMO naming convention variants

OUTPUT

Same text files but with the headers renamed in both the input matrix
and observations/colmeta files.

If you only pass the --input_directory it assumes the files are called 'expression.tab',
'observations.tab' and 'genes.tab'

{input_directory}/expression.tab
{input_directory}/observations.tab
{input_directory}/genes.tab

If you also pass --file_prefix we're assuming NeMO file naming conventions and it will
instead look for:

{input_directory}/{file_prefix}_DataMTX.tab
{input_directory}/{file_prefix}_COLmeta.tab
{input_directory}/{file_prefix}_ROWmeta.tab
"""

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import os
import shutil
import sys


def main():
    parser = argparse.ArgumentParser(description='Replaces numeric headers by prepending values to make them strings')

    parser.add_argument('-i', '--input_directory', type=str, required=True, help='Path to an input directory containing the three tab files' )
    parser.add_argument('-f', '--file_prefix', type=str, required=False, help='Prefix portion for each of the three input files' )
    parser.add_argument('-c', '--column_prefix', type=str, required=False, default='X', help='Text which will get pre-pended to each column header' )
    args = parser.parse_args()

    if args.file_prefix:
        mtx_file = "{0}/{1}_DataMTX.tab".format(args.input_directory, args.file_prefix)
        obs_file = "{0}/{1}_COLmeta.tab".format(args.input_directory, args.file_prefix)
        var_file = "{0}/{1}_ROWmeta.tab".format(args.input_directory, args.file_prefix)
    else:
        mtx_file = args.input_directory + '/expression.tab'
        obs_file = args.input_directory + '/observations.tab'
        var_file = args.input_directory + '/genes.tab'

    for filename in [mtx_file, obs_file, var_file]:
        if os.path.exists("{0}.new".format(filename)):
            print("ERROR: this script write temporary '.new' file names while working but one already exists. Quitting so it isn't stomped", file=sys.stderr)
            sys.exit(1)

    # key:old value,  value: new value
    replaced_columns = dict()
    mtx_new_filename = "{0}.new".format(mtx_file)
    obs_new_filename = "{0}.new".format(obs_file)
        
    with open(mtx_file) as fh:
        line_count = 0
        out_fh = open(mtx_new_filename, 'wt')
        for line in fh:
            if line_count:
                out_fh.write(line)
            else:
                line = line.rstrip()
                cols = line.split("\t")
                if len(cols) > 1:
                    out_fh.write("{0}\t".format(cols[0]))
                    new_cols = list()
                    for col in cols[1:]:
                        new_cols.append("{0}{1}\t".format(args.column_prefix, col))
                        replaced_columns[str(col)] = "{0}{1}".format(args.column_prefix, col)
                    out_fh.write("{0}\n".format("\t".join(new_cols)))
                    out_fh.write("\n")
                else:
                    out_fh.write(line)

            line_count += 1
        out_fh.close()

    with open(obs_file) as fh:
        line_count = 0
        out_fh = open(obs_new_filename, 'wt')
        for line in fh:
            if line_count:
                line = line.rstrip()
                cols = line.split("\t")
                if cols[0] in replaced_columns:
                    out_fh.write("{0}\t".format(replaced_columns[cols[0]]))
                    out_fh.write("{0}\n".format("\t".join(cols[1:])))
                else:
                    raise Exception("ERROR: column was found in the obs index but not the matrix headers")
            else:
                out_fh.write(line)

            line_count += 1

        out_fh.close()

    ## If we got this far, move the temp files over the originals.
    shutil.move(mtx_new_filename, mtx_file)
    shutil.move(obs_new_filename, obs_file)
        
        
if __name__ == '__main__':
    main()







