#!/usr/bin/env python3

"""
Adds a gene symbol column to a tab delimited text file.

Assumes ensembl ID is the first column.
"""

import argparse
import mysql.connector

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created.')
    parser.add_argument('-org', '--organism_id', type=int, required=True, help='Organism ID')
    parser.add_argument('-r', '--release_num', type=int, required=True, help='Ensembl release number')
    parser.add_argument('-u', '--user', type=str, required=True, help='Database user.')
    parser.add_argument('-p', '--password', type=str, required=True, help='Database password.')
    args = parser.parse_args()

    # connect to db
    cnx = mysql.connector.connect(user=args.user, password=args.password, host='localhost', database='gear_portal')
    cursor = cnx.cursor()

    # get organism to filter out gene table
    qry = "SELECT ensembl_id, gene_symbol FROM gene WHERE organism_id = %s AND ensembl_release = %s"
    cursor.execute(qry, (args.organism_id, args.release_num))

    gene_syms = dict()
    for (id, gs) in cursor:
        gene_syms[id] = gs

    cursor.close()
    cnx.close()

    ofh = open(args.output_file, 'wt')

    missing_count = 0
    
    line_num = 0
    for line in open(args.input_file):
        line = line.rstrip()
        line_num += 1

        if line_num == 1:
            print("{0}\tgene_symbol".format(line), file=ofh)
        else:
            ensembl_id = line.split("\t")[0]
            if ensembl_id in gene_syms:
                print("{0}\t{1}".format(line, gene_syms[ensembl_id]), file=ofh)
            else:
                missing_count += 1
                print("{0}\t".format(line), file=ofh)

    print("Gene symbols added.  There were {0} ensembl_ids in the input which weren't found for mapping".format(missing_count))
    
if __name__ == '__main__':
    main()
