#!/usr/bin/env python3

"""

The observations.tab file contains entries like:

observations    cell_type       replicate       seurat_clusters
D3-1_AAATGGAGTAGCTCGC-1 mESCs   1       1
D3-1_AACCAACCATGTTCAG-1 mESCs   1       1
D3-1_AAGTGAACACTTCTCG-1 transitional    1       3
D3-1_ACAACCATCCGTAATG-1 ectoderm        1       4
...

And expression.tab starts with:

D3-1_AAATGGAGTAGCTCGC-1 D3-1_AACCAACCATGTTCAG-1 D3-1_AAGTGAACACTTCTCG-1 ...

The rows in observations.tab should correspond to the columns in expression.tab, and this
script checks that this is true.

"""

import argparse
import math
import os


def main():
    parser = argparse.ArgumentParser( description='Validates the correlation between expression.tab and observations.tab files.')
    parser.add_argument('-e', '--expression_file', type=str, required=True, help='Path to the expression.tab file' )
    parser.add_argument('-o', '--observations_file', type=str, required=True, help='Path to the observations.tab file' )
    args = parser.parse_args()

    line_num = 0
    obs = list()
    
    ## first read in the observations
    for line in open(args.observations_file):
        line = line.rstrip()
        line_num += 1

        if line_num > 1:
            cols = line.split("\t")
            obs.append(cols[0])

    print("Observations count: {0}".format(len(obs)))

    line_num = 0
    expr = None
    
    for line in open(args.expression_file):
        line = line.rstrip()
        line_num += 1

        if line_num == 1:
            expr = line.split("\t")

    print("Expression column count: {0}".format(len(expr)))

    ## Now actually compare
    i = 0
    for oval in obs:
        if oval != expr[i]:
            raise Exception("obs value {0} didn't match expr value {1} at index {2} (0-based)".format(oval, expr[i], i))
        
        i += 1
        

    print("obs and expr values match properly")

if __name__ == '__main__':
    main()







