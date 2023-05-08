#!/usr/bin/env python3

"""
Print best match ensembl release for set of ensembl ids from file, count, and organism id.

Assumes Ensembl ID or gene symbol is the first column of the input file
"""

import argparse
import os, sys

lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'lib')
sys.path.append(lib_path)
import geardb

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read.')
    parser.add_argument('-org', '--organism_id', type=int, required=True, help='Organism ID')
    parser.add_argument("--silent", action="store_true", required=False, help="Suppress printing. Only return best release version")
    args = parser.parse_args()
    find_best_ensembl_release_match(args.input_file, args.organism_id, args.silent)

def find_best_ensembl_release_match(input_file, organism_id, silent=False):
    qry_ids = list()
    id_ref = dict()

    line_num = 0
    for line in open(input_file):
        line = line.rstrip()
        line_num += 1

        if line_num > 1:
            cols = line.split("\t")
            qry_ids.append(cols[0])

    # connect to db
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    key_type = get_key_type(input_file)

    # get organism to filter out gene table
    qry = "SELECT {0}, ensembl_release FROM gene WHERE organism_id = %s".format(key_type)
    cursor.execute(qry, (organism_id,))

    for (id, release_num) in cursor:
        if release_num not in id_ref:
            id_ref[release_num] = dict()

        id_ref[release_num][id] = 1

    best_release = None
    best_count = 0
    best_missing_count = 0

    for release_num in sorted(id_ref):
        cov_count, empty_count = get_count(qry_ids, id_ref[release_num])
        if not silent:
            print("release:{0}\tcount:{1}".format(release_num, cov_count))

        if cov_count > best_count:
            best_count = cov_count
            best_release = release_num
            best_missing_count = empty_count
    if not silent:
        print("The best release is {0} with {1} of {2} genes unaccounted for.  Of these, {3} were empty".format(
            best_release, len(qry_ids) - cov_count, len(qry_ids), best_missing_count))

    cursor.close()
    cnx.close()

    return best_release

def get_count(queries, refs):
    c = 0
    m = 0

    for q in queries:
        if q in refs:
            c += 1

        elif q is None or q == '':
            m += 1

    return c, m

def get_key_type(ifile):
    """
    We want to magically allow for the first column to be either gene symbols or Ensembl gene IDs.
    This determines the type by checking the first column of the first 10 rows (skipping the header)
    and returns 'ensembl_id' if all start with ENS, else 'gene_symbol'
    """
    line_num = 0
    for line in open(ifile):
        line_num += 1
        if line_num == 1: continue
        if line_num > 10: break

        key = line.split("\t")[0]

        if not key.startswith('ENS'):
            return 'gene_symbol'

    return 'ensembl_id'

if __name__ == '__main__':
    main()









