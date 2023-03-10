#!/opt/bin/python3

"""
Given an h5ad file with only gene symbols, find best matching ensembl release, then
collect ensembl ids for the corresponding gene symbol.

This version is a little different because it relies on file manipulation and the following directories
to exist in the cwd for it to work:

$ mkdir mapped unmapped merged

Kept failing to maintain the var dataframe when using pandas/scanpy, so this method was created.

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
import shutil

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
    parser.add_argument('-uec', '--use_ensembl_column', type=str,
                        required=False, help="If the datafile already has an Ensembl column it wasn't indexed with, pass its name here.")
    

    args = parser.parse_args()

    adata = sc.read(args.input_file)
    (_, n_genes) = adata.shape

    ensembl_releases = [84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94]
    #ensembl_releases = [84, 85]
    
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

    adata = adata[:, ~duplicated_genes]

    print("\nOriginal loaded adata\n")
    print(adata)

    # If an ensembl column was provided all we have to do is reindex on that and rewrite the file
    if args.use_ensembl_column:
        print("INFO: re-indexing file on column ({0})".format(args.use_ensembl_column))
        adata.var = adata.var.set_index(args.use_ensembl_column)
        
        if not args.read_only:
            adata.write(args.output_file)

        # And we're done
        print('############################################################')
        return True

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

    # What if we didn't find any? Throw error and fail
    if best_release is None:
        raise Exception("ERROR: No mappings found (check that the proper organism was passed?)")
    
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
        var=ensembl_id_var)

    ## Now combine the unmapped dataframe with this one, first making the needed edits
    if 'gene_symbol' in adata_not_present.var.columns:
        adata_not_present.var = adata_not_present.var.rename(columns={"gene_symbol": "gene_symbol_original"})

    # Splitting code over multiple lines requires a "\" at the end.
    adata_unmapped_var = adata_not_present.var.reset_index() \
        .rename(columns={ \
                    #adata_not_present.var.index: "ensembl_id",
                    "genes": "gene_symbol" \
                    }) \
        .set_index(args.id_prefix + adata_not_present.var.index.astype(str))

    adata_unmapped = ad.AnnData(
        X=adata_not_present.X,
        obs=adata_not_present.obs,
        var=adata_unmapped_var
    )
    adata_unmapped.var.index.name = "ensembl_id"

    print("ADATA UNMAPPED.VAR")
    print(adata_unmapped.var)

    ## write the mapped to a set of files
    adata_with_ensembl_ids.transpose().write_csvs('mapped', sep="\t", skip_data=False)
    adata_unmapped.transpose().write_csvs('unmapped', sep="\t", skip_data=False)

    shutil.copyfile('mapped/obs.csv', 'merged/obs.tab')

    # place the var and X and then add to each
    shutil.copyfile('mapped/var.csv', 'merged/var.tab')
    shutil.copyfile('mapped/X.csv', 'merged/X.tab')

    # place the mapped file and then
    append_file('unmapped/obs.csv', 'merged/obs.tab', True)
    append_file('unmapped/X.csv', 'merged/X.tab', False)

    # parse the text files for a new adata (the transposition switches what we'd call var and obs)
    adata = sc.read('merged/X.tab', sep='\t', cache=False).transpose()
    obs = pd.read_table('merged/var.tab', sep='\t', index_col=0, header=0)
    var = pd.read_table('merged/obs.tab', sep='\t', index_col=0, header=0)

    for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
        if str_type in obs.columns:
            obs[str_type] = pd.Categorical(obs[str_type])

    for num_type in ['replicate', 'time_point_order']:
        if num_type in obs.columns:
            obs[num_type] = pd.to_numeric(obs[num_type])

    # Assign genes and observations to AnnData object
    adata.var = var
    adata.obs = obs

    if not args.read_only:
        adata.write(args.output_file)
    print('############################################################')


def append_file(src, dest, skip_header):
    ofh = open(dest, 'a')

    line_num = 0

    for line in open(src):
        if skip_header:
            if line_num > 0:
                ofh.write(line)

            line_num += 1
        else:
            ofh.write(line)

    ofh.close()

if __name__ == '__main__':
    main()
