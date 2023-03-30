#!/usr/bin/env python3

"""
"""

import argparse
import configparser
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import os
import sys
import mysql.connector
from mysql.connector import errorcode

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-org', '--organism', type=int, required=True, help='Organism ID to use')
    parser.add_argument('-er', '--ensembl_release', type=int, required=True, help='Ensembl release ID to use')
    args = parser.parse_args()

    ensembl_genes_df = get_df_index(args)
    adata = sc.read_h5ad(args.input_file)

    # Some genes are duplicated in the index,
    # we make them unique here, however we lose
    # these genes when merging with our ensembl ids
    # because the new unique names do not match
    # genes from the sql query.
    adata.var_names_make_unique()

    # now rename the index to that column
    adata.var.index.rename('gene_symbol', inplace=True)

    # Do an inner join between adata's current var
    # and the ensembl_id/gene_symbol table fetched
    # from gene table
    merged_var = (
        adata
            .var
            .reset_index()
            .merge(ensembl_genes_df, on='gene_symbol')
            .set_index('ensembl_id')
    )

    # Filter adata based on gene's that exist in our sql query
    gene_filter = adata.var.index.isin(merged_var.gene_symbol)
    adata = adata[:, gene_filter]

    # Currently creating a new AnnData object because of
    # trouble getting adata.var = merged_var to persist
    adata = ad.AnnData(adata.X, obs=adata.obs, var=merged_var)

    adata.write(args.output_file)

def get_df_index(args):
    """
    Returns Pandas dataframe with Ensembl IDs and gene symbols.
    """
    config = configparser.ConfigParser()
    config.read('../gear.ini')

    cnx = None
    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'])
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password", file=sys.stderr)
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist", file=sys.stderr)
        else:
            print(err, file=sys.stderr)

    cursor = cnx.cursor()

    query = """
        SELECT ensembl_id, gene_symbol
          FROM gene
         WHERE organism_id = %s
           AND ensembl_release = %s
    """
    cursor.execute(query, (args.organism, args.ensembl_release))

    df = pd.DataFrame(cursor.fetchall(), columns=cursor.column_names)
    # Some Gene's are mapped to multiple ensembl IDs, so we
    # drop the duplicate's and take the first Ensembl ID for now.
    df.drop_duplicates('gene_symbol', inplace=True)
    df.set_index('gene_symbol', inplace=True)

    cursor.close()
    cnx.close()

    return df

if __name__ == '__main__':
    main()
