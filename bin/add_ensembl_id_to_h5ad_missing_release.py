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

    print("\nOriginal loaded adata\n")
    print(adata)

    # If adata.var is an empty dataframe, make note that the index is the original gene symbol column
    # Ensures the `adata_unmapped_var` rename aligns with the original gene symbol column in adata.var
    orig_gene_column = "genes"
    if adata.var.empty:
        orig_gene_column = "index"

    for release in ensembl_releases:
        print("INFO: comparing with ensembl release: {0} ... ".format(release), end='')
        cursor.execute(query, (args.organism, release))

        df = pd.DataFrame(cursor.fetchall(), columns=cursor.column_names)

        # Query can return different ensembl ids for a gene symbol,
        # we want to drop duplicates so we have only one ensembl id
        # per gene symbol
        df = df.drop_duplicates(subset=['gene_symbol'])
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
    adata_not_present = adata[:, ~genes_present_filter]

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

    print("ENSEMBL_ID_VAR")
    print(ensembl_id_var)

    # Currently creating a new AnnData object because of
    # trouble getting adata.var = merged_var to persist
    adata_with_ensembl_ids = ad.AnnData(
        adata_present.X,
        obs=adata_present.obs,
        var=ensembl_id_var,
        # May not use these directly, but need to pass them through to preserve them
        # Note: there are other fields, like layers, that could be added here if needed
        obsm=adata_present.obsm,
        obsp=adata_present.obsp,
        varm=adata_present.varm,
        varp=adata_present.varp,
        uns=adata_present.uns
        )

    ## Now combine the unmapped dataframe with this one, first making the needed edits
    if 'gene_symbol' in adata_not_present.var.columns:
        adata_not_present.var = adata_not_present.var.rename(columns={"gene_symbol": "gene_symbol_original"})

    # Splitting code over multiple lines requires a "\" at the end.
    adata_unmapped_var = adata_not_present.var.reset_index() \
        .rename(columns={ \
                    #adata_not_present.var.index: "ensembl_id",
                    # In order to concatenate the unmapped var with the mapped var, we need to ensure the
                    # column names are not the same as the mapped var. The "concat" function treats columns
                    # as a set of values, so for a strategy like "unique" which requires all values
                    # to be the same, it would see two different sets of values and not merge them.
                    orig_gene_column: "gene_symbol_unmapped" \
                    }) \
        .set_index(args.id_prefix + adata_not_present.var.index.astype(str))

    adata_unmapped = ad.AnnData(
        X=adata_not_present.X,
        obs=adata_not_present.obs,
        var=adata_unmapped_var,
        obsm=adata_not_present.obsm,
        obsp=adata_not_present.obsp,
        varm=adata_not_present.varm,
        varp=adata_not_present.varp,
        uns=adata_not_present.uns
    )
    adata_unmapped.var.index.name = "ensembl_id"

    print("ADATA UNMAPPED.VAR")
    print(adata_unmapped.var)

    # Concatenate the two AnnData objects, which unfortunately duplicates the observations
    adata = ad.concat([adata_with_ensembl_ids, adata_unmapped], join="outer", merge="unique", uns_merge="unique", label="dataset")

    #print("ADATA CONCAT")
    #print(adata)
    #print("VAR CONCAT")
    #print(adata.var)
    #print("OBS_CONCAT")
    #print(adata.obs)

    # Backfill the gene symbol column with the unmapped gene symbols
    adata.var['gene_symbol'] = adata.var['gene_symbol'].combine_first(adata.var['gene_symbol_unmapped'])

    # Drop the unmapped gene symbol column
    adata.var = adata.var.drop(columns=['gene_symbol_unmapped'])

    # There should be duplicated indexes in the adata.obs dataframe, so drop the duplicates
    # We could alternatively use `duplicates()` but this feels safer
    adata = adata[adata.obs.dataset == "0"]
    adata.obs = adata.obs.drop(columns=['dataset'])



    print("FINAL ADATA")
    print(adata)
    print('VAR\n')
    print(adata.var)
    print("OBS\n")
    print(adata.obs)


    if not args.read_only:
        adata.write(args.output_file)
    print('############################################################')


if __name__ == '__main__':
    main()
