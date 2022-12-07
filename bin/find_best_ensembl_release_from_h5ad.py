#!/usr/bin/env python3

"""
Print best match ensembl release for set of ensembl ids from file, count, and organism id.

Assumes Ensembl ID is the first column of the input file
"""

import argparse
import sys
import os
import scanpy as sc
import mysql.connector

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)

import geardb

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read.')
    parser.add_argument('-org', '--organism_id', type=int, required=True, help='Organism ID')
    args = parser.parse_args()

    adata = sc.read(args.input_file)
    qry_ids = adata.var.index.tolist()
    id_ref = dict()

    first_index = adata.var.first_valid_index()
    key_type = get_key_type(first_index)
    
    # connect to db
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    # get organism to filter out gene table
    if key_type == 'ensembl_id':
        qry = "SELECT ensembl_id, ensembl_release FROM gene WHERE organism_id = %s"
    elif key_type == 'gene_symbol':
        qry = "SELECT gene_symbol, ensembl_release FROM gene WHERE organism_id = %s"
    else:
        raise Error("Couldn't guess key type based on first key ({0})".format(first_index))

    cursor.execute(qry, (args.organism_id,))

    for (id, release_num) in cursor:
        if release_num not in id_ref:
            id_ref[release_num] = dict()

        id_ref[release_num][id] = 1;

    best_release = None
    best_count = 0

    for release_num in id_ref:
        cov_count = get_count(qry_ids, id_ref[release_num])
        print("release:{0}\tcount:{1}".format(release_num, cov_count))

        if cov_count > best_count:
            best_count = cov_count
            best_release = release_num

    print("The best release is {0} with {1} of {2} genes unaccounted for.".format(
        best_release, len(qry_ids) - cov_count, len(qry_ids)))

    cursor.close()
    cnx.close()

def get_count(queries, refs):
    c = 0

    for q in queries:
        if q in refs:
            c += 1

    return c

def get_key_type(id):
    """
    We want to magically allow for the first column to be either gene symbols or Ensembl gene IDs.

    This checks to see if the identifier begins with ENS and returns either 'gene_symbol' 
    or 'ensembl_id'
    """
    if not id.startswith('ENS'):
        return 'gene_symbol'

    return 'ensembl_id'

if __name__ == '__main__':
    main()
