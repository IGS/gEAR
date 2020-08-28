#!/usr/bin/env python3

"""
Our supported 3-tab Excel sheet format has one tab called 'genes' which can have a variety of
columns.  The ones we care about here are

   id - Ensembl IDs
   gene_symbol - Individual gene symbols for each gene

This script checks that the gene symbol listed in the column matches the reference annotation
we have stored in the database.
"""

import argparse
import os
import pandas as pd
from pandas import ExcelFile
import re
import sys

lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)

import gear.db

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument("-v", "--versioned", help="Input ENSEMBL IDs are versioned (with .N extensions)", action="store_true")
    args = parser.parse_args()

    #connect to gEAR MySQL
    mysql_cnx = gear.db.MySQLDB().connect()
    cursor = mysql_cnx.cursor()
    ensembl_idx = get_gene_symbol_lookup(cursor=cursor, use_versions=args.versioned)

    #for id in ensembl_idx:
    #    print("{0}\t{1}".format(id, ensembl_idx[id]))

    genes_df = pd.read_excel(args.input_file, sheet_name='genes', index_col='id', converters={'gene_symbol': str})

    for i in genes_df.index:
        gs = genes_df['gene_symbol'][i]
        id = i

        m = re.match('^\d+$', gs)
        if m:
            #print("This gene symbol appeared to be only digits: {0}".format(gs))

            if id in ensembl_idx:
                #print("\tLooks like I could replace it with any of: {0}".format(ensembl_idx[id]))
                print("{0}\t{1}".format(id, ensembl_idx[id][0]))

    #for gs in genes_df['gene_symbol']:
    #    m = re.match('^\d+$', gs)
    #    if m:
    #        print("This gene symbol appeared to be only digits: {0}".format(gs))

    cursor.close()


def get_gene_symbol_lookup(cursor=None, use_versions=False):
    # Returns a dictionary where key=ensembl_id and value is gene symbol
    idx = dict()

    if use_versions == True:
        qry = "SELECT ensembl_version, gene_symbol FROM gene ORDER BY ensembl_release DESC"
    else:
        qry = "SELECT ensembl_id, gene_symbol FROM gene ORDER BY ensembl_release DESC"
    
    cursor.execute(qry)

    for [ensembl_id, gene_symbol] in cursor:
        if ensembl_id in idx:
            if gene_symbol not in idx[ensembl_id]:
                idx[ensembl_id].append(gene_symbol)
        else:
            idx[ensembl_id] = [gene_symbol]

    return idx


if __name__ == '__main__':
    main()
