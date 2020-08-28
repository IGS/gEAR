#!/usr/bin/env python3

"""
ALTER TABLE gene ADD column molecule VARCHAR(100) after organism_id;
ALTER TABLE gene ADD column start INT after molecule;
ALTER TABLE gene ADD column stop INT after start;


Author: Joshua Orvis (jorvis AT gmail)
"""

import argparse
import os
import sys
import re
import gzip
import mysql.connector
import configparser
from biocode import utils

import loaderutils

import gear

from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=False, help='Path to an input GBK file' )
    parser.add_argument('-l', '--input_list', type=str, required=False, help='Path to an input GBK list file' )
    parser.add_argument('-oid', '--organism_id', type=str, required=True, help='Organism ID being loaded' )

    files = list()
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read('gear.ini')

    if args.input_file is not None:
        files.append(args.input_file)
    elif args.input_list is not None:
        files = utils.read_list_file(args.input_list)
    else:
        raise Exception("ERROR: You must pass either -i or -l options")

    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'])
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)

    cursor = cnx.cursor()
    genes_by_ensembl_id = gear.cache_genes_by_ensembl_id(cursor)
    genes_by_sym = loaderutils.cache_genes_by_primary_gene_symbol(cursor, lower=True)

    for file in files:
        recompress = False
        
        # decompress the file if needed.
        if file.endswith(".gz"):
            os.system("gunzip {0}".format(file))
            recompress = True
            m = re.match("(.+\.([A-Za-z0-9]+)\.dat)\.gz", file)
            if m:
                file = m.group(1)
                chromosome = m.group(2)
            else:
                raise Exception("This should not have happened.  Failed to regex file path")
        else:
            recompress = False
        
        # each gb_record is a SeqRecord object
        print("Processing file {0} ...".format(file))
        for gb_record in SeqIO.parse(open(file, "r"), "genbank"):
            mol_id = gb_record.name

            # each feat is a SeqFeature object
            for feat in gb_record.features:
                if 'gene' in feat.qualifiers:
                    ensembl_id = feat.qualifiers['gene'][0]

                    # remove any versioning
                    if '.' in ensembl_id:
                        ensembl_id = ensembl_id.split('.')[0]

                    if ensembl_id not in genes_by_ensembl_id:
                        continue
                    
                if feat.type == 'gene':
                    gene_id = genes_by_ensembl_id[ensembl_id]
                    
                    start = feat.location.start.position
                    stop  = feat.location.end.position

                    print("{0}\t{1}\t{2}".format(ensembl_id, start, stop))

        # now recompress the file, if necessary
        if recompress:
            os.system("gzip {0}".format(file))

    cnx.commit()
    cursor.close()
    cnx.close()


if __name__ == '__main__':
    main()







