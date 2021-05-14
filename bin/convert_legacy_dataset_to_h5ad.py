#!/opt/bin/python3

"""

This script takes a dataset tab file created by the legacy gEAR uploader and converts
it to an AnnData H5AD file.  Previously, the column header format convention provided
some metadata for the columns.  These conventions and how they are transformed are:

Bar and line plots:
   Condition--Replicate_[pval|sd]
   Untreated--P_0
   Treated--P_1_sd
   Treated--P_1_pval

Violin plot
   Group--Label
   OHC--Cell1
   OHC--Cell2

SVG
   Label

"""

import argparse
import gzip
import os
import re
import sys

import mysql.connector

import anndata
import pandas as pd
import scanpy.api as sc

sys.path.append("{0}/../lib".format(os.path.dirname(sys.argv[0])))
import geardb

def main():
    parser = argparse.ArgumentParser( description='Annotation loader for the gEAR')
    parser.add_argument('-i', '--input_tab', type=str, required=True, help='Path to an input file to be read')
    parser.add_argument('-o', '--output_h5ad', type=str, required=True, help='Path to the output file to be created')
    
    args = parser.parse_args()

    cnx = geardb.Connection()
    cursor = geardb.Connection().get_cursor()
    
    cursor.close()
    cnx.close()
            
if __name__ == '__main__':
    main()







