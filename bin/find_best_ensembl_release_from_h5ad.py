#!/usr/bin/env python3

"""
Print best match ensembl release for set of ensembl ids from file, count, and organism id.
Assumes Ensembl ID is the first column of the input file
"""

import argparse
import scanpy as sc
import mysql.connector
import geardb

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read.')
    parser.add_argument('-org', '--organism_id', type=int, required=True, help='Organism ID')
    args = parser.parse_args()

    adata = sc.read(args.input_file)
    qry_ids = adata.var.index.tolist()
    id_ref = dict()

    # connect to db
    #cnx = mysql.connector.connect(user=args.user, password=args.password, host='localhost', database='gear_portal')
    cnx = geardb.Connection()
    #cursor = cnx.cursor()
    cursor = cnx.get_cursor()

    # get organism to filter out gene table
    qry = "SELECT ensembl_id, ensembl_release FROM gene WHERE organism_id = %s"
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

if __name__ == '__main__':
    main()
