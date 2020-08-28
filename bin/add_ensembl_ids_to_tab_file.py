#!/opt/bin/python3

"""

"""

import argparse
import cgi, json
import mysql.connector
import os, sys

lib_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 'lib')
sys.path.append(lib_path)
import geardb


def main():
    parser = argparse.ArgumentParser( description='See the script name?  It does that.')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-org', '--organism', type=int, required=True, help='Organism ID to use')
    parser.add_argument('-er', '--ensembl_release', type=int, required=True, help='Ensembl release ID to use')
    args = parser.parse_args()
    
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    query = """
        SELECT ensembl_id, gene_symbol
          FROM gene
         WHERE organism_id = %s
           AND ensembl_release = %s
    """

    gs_idx = dict()
    duplicate_gs = 0
    gs_not_found = 0
    
    cursor.execute(query, (args.organism, args.ensembl_release))
    for row in cursor:
        ensembl_id = row[0]
        gs = row[1]

        if gs not in gs_idx:
            gs_idx[gs] = ensembl_id
        else:
            duplicate_gs += 1

    print("INFO: There were {0} instances where a gene symbol was associated with more than one Ensembl ID".format(duplicate_gs), file=sys.stderr)

    cursor.close()
    cnx.close()
            
    line_num = 0

    ofh = open(args.output_file, 'w')
    
    for line in open(args.input_file):
        line_num += 1
        line = line.rstrip()
        
        if line_num == 1:
            print("Ensembl_ID\t" + line, file=ofh)
            continue

        cols = line.split("\t")
        
        if cols[0] in gs_idx:
            cols.insert(0, gs_idx[cols[0]])
            print("\t".join(cols), file=ofh)
        else:
            print("Gene not found: ({0})".format(cols[0]))
            gs_not_found += 1

    print("INFO: There were {0} instances of rows skipped in output because their gene symbol wasn't found in the index".format(gs_not_found), file=sys.stderr)

if __name__ == '__main__':
    main()
