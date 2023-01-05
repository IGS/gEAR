#!/opt/bin/python3

"""
Given an h5ad file with only gene symbols, find best matching ensembl release, then
collect ensembl ids for the corresponding gene symbol. Some information loss is possible
with the duplication of gene symbols.
"""

import argparse
import pandas as pd
import anndata as ad
import scanpy as sc
import cgi
import json
import mysql.connector
import sys
import os

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)

import geardb

def main():
    parser = argparse.ArgumentParser(
        description='See the script name?  It does that.')
    parser.add_argument('-i', '--input_file', type=str,
                        required=True, help='Path to an h5ad file to be read.')
    parser.add_argument('-o', '--output_file', type=str,
                        required=True, help='Output file to be written.')
    parser.add_argument('-org', '--organism', type=int,
                        required=True, help='Organism ID.')
    parser.add_argument('-r', '--read_only', type=bool, default=False,
                        help='If true, only print and do not write to output.')

    args = parser.parse_args()

    adata = sc.read(args.input_file)
    (_, n_genes) = adata.shape

    ensembl_releases = [84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94]

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    query = """
        SELECT ensembl_id, gene_symbol
          FROM gene
         WHERE organism_id = %s
           AND ensembl_release = %s
    """

    best_release = None
    best_count = 0
    best_df = None

    # There are some cases where there are duplicate gene symbol,
    # we do not want to keep these so we drop them here in order
    # preserve obs shape when joining ensembl ids and resassigning
    # var later on.
    duplicated_genes = adata.var.index.duplicated()
    print('############################################################')
    print(f"File: {args.input_file}")
    print(f"Duplicated Genes: {duplicated_genes.sum()}")

    # print(f"The file, {args.input_file} has {duplicated_genes.sum()} duplicate genes. These will be dropped.")
    adata = adata[:, ~duplicated_genes]
    var = adata.var

    for release in ensembl_releases:
        print("INFO: comparing with ensembl release: {0} ... ".format(release), end='')
        cursor.execute(query, (args.organism, release))

        df = pd.DataFrame(cursor.fetchall(), columns=cursor.column_names)

        # Query can return different ensembl ids for a gene symbol,
        # we want to drop duplicates so we have only one ensembl id
        # per gene symbol
        df = df.drop_duplicates('gene_symbol')
        df = df.set_index('gene_symbol')

        merged_df = var.join(df, how='inner')
        (row_count, _) = merged_df.shape

        print(" found {0} matches".format(row_count))

        if row_count > best_count:
            best_count = row_count
            best_release = release
            best_df = merged_df

    print(f"Best release: {best_release}")
    print(f"Matches for release: {best_count}")
    print(f"Original # Genes: {n_genes}")
    print(f"Genes lost: {n_genes - best_count}")
    # Now we have our best release and ensembl ids for those gene symbols,
    # we want to make sure we filter our anndata object to reflect only the
    # gene symbols we have ids for
    gene_filter = var.index.isin(best_df.index)
    adata = adata[:, gene_filter]

    # If the data already had a 'gene symbol' let's rename it
    if 'gene_symbol' in best_df.columns:
        print("WARN: Found gene_symbol column already in dataset, renaming to gene_symbol_original")
        best_df = best_df.rename(columns={"gene_symbol": "gene_symbol_original"})

    ensembl_id_var = (
        best_df
        .reset_index()
        .rename(columns={
            "index": "gene_symbol"
        })
        .set_index('ensembl_id')
    )

    # Currently creating a new AnnData object because of
    # trouble getting adata.var = merged_var to persist
    adata_with_ensembl_ids = ad.AnnData(
        adata.X,
        obs=adata.obs,
        var=ensembl_id_var)
    print('VAR\n')
    print(adata_with_ensembl_ids.var.head())
    print("OBS\n")
    print(adata_with_ensembl_ids.obs.head())
    if not args.read_only:
        adata_with_ensembl_ids.write(args.output_file)
    print('############################################################')


if __name__ == '__main__':
    main()
