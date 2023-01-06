#!/opt/bin/python3

"""
Given an h5ad file with only gene symbols, find best matching ensembl release, then
collect ensembl ids for the corresponding gene symbol. 

Fake ensembl IDs will be created for those genes whose gene symbols don't map to
a known ensembl ID.
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
    parser.add_argument('-idp', '--id_prefix', type=str,
                        required=True, help="Fake Ensembl IDs will be generated for any rows whose gene symbols don't match, using this prefix.")

    args = parser.parse_args()

    adata = sc.read(args.input_file)
    (_, n_genes) = adata.shape

    #ensembl_releases = [84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94]
    ensembl_releases = [84, 85]

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

    print("\nOriginal loaded adata\n")
    print(adata)

    for release in ensembl_releases:
        print("INFO: comparing with ensembl release: {0} ... ".format(release), end='')
        cursor.execute(query, (args.organism, release))

        df = pd.DataFrame(cursor.fetchall(), columns=cursor.column_names)

        # Query can return different ensembl ids for a gene symbol,
        # we want to drop duplicates so we have only one ensembl id
        # per gene symbol
        df = df.drop_duplicates('gene_symbol')
        df = df.set_index('gene_symbol')

        merged_df = adata.var.join(df, how='inner')
        (row_count, _) = merged_df.shape

        print(" found {0} matches".format(row_count))

        if row_count > best_count:
            best_count = row_count
            best_release = release
            best_df = merged_df

    print(f"\nBest release: {best_release}")
    print(f"Matches for release: {best_count}")
    print(f"Original # Genes: {n_genes}")
    print(f"Genes lost: {n_genes - best_count}\n")
    
    # Now we have our best release and ensembl ids for those gene symbols,

    # Get separate adata for those where the gene symbols were mapped and where they weren't
    genes_present_filter = adata.var.index.isin(best_df.index)
    adata_present = adata[:, genes_present_filter]
    adata_unmapped = adata[:, ~genes_present_filter]

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
        adata_present.X,
        obs=adata_present.obs,
        var=ensembl_id_var)

    print('\nVAR with ensembl ids\n')
    print(adata_with_ensembl_ids.var.head())

    ## Now combine the unmapped dataframe with this one, first making the needed edits
    if 'gene_symbol' in adata_unmapped.var.columns:
        adata_unmapped.var = adata_unmapped.var.rename(columns={"gene_symbol": "gene_symbol_original"})

    print("]nVAR of adata_unmapped\n")
    print(adata_unmapped.var_names)

    adata_unmapped_var = adata_unmapped.var.reset_index()
    adata_unmapped_var = adata_unmapped_var.set_index(args.id_prefix + adata_unmapped.var.index.astype(str))
    adata_unmapped_var = adata_unmapped_var.rename(columns={
        adata_unmapped_var.columns[0]: "ensembl_id",
        "genes": "gene_symbol"
    })

    print("]nVAR of adata_unmapped_var\n")
    print(adata_unmapped_var.head(5))
    
    adata_unmapped = ad.AnnData(
        X=adata_unmapped.X,
        obs=adata_unmapped.obs,
        var=adata_unmapped_var
    )

    adata = adata_with_ensembl_ids.concatenate(adata_unmapped2)
    
    print('VAR\n')
    print(adata.var.head())
    print("OBS\n")
    print(adata.obs.head())
    if not args.read_only:
        adata.write(args.output_file)
    print('############################################################')


if __name__ == '__main__':
    main()
