#!/usr/bin/env python3

"""
load_komp_repository_data.py - Load gene synonyms and external links scraped from the KOMP repository.

Reads komp_repository_data.p (a pickle written by fetch_komp_repository_data.py) and, for each
gene symbol already in the gene table, adds its KOMP synonyms to gene_symbol as non-primary
symbols and its MGI, IGTC, IMSR, BioGPS and GEO links to gene_urls. The gene_urls table is not
defined in create_schema.sql, and nothing on the website reads it. Reads database settings from
../gear.ini, so run it from bin/.
"""

import argparse
import mysql.connector
import configparser
import pickle
import os
import re
import sys
sys.path.insert(0, '../')
import loaderutils


def main():
    config = configparser.ConfigParser()
    config.read('../gear.ini')

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

    gene_cache = loaderutils.cache_genes_by_primary_gene_symbol(cursor)

    KOMP_BASE_URL = 'https://www.komp.org'
    PICKLE_FILE = 'komp_repository_data.p'

    dstore = pickle.load(open(PICKLE_FILE, 'rb'))

    add_gs = "INSERT INTO gene_symbol (gene_id, label, is_primary) VALUES (%s, %s, 0)"

    for gene_sym in dstore:
        if gene_sym in gene_cache:
            gene_id = gene_cache[gene_sym]
            print("GENE SYM: {0}".format(gene_sym))
            print("\tSYNONYMS:")

            for syn in dstore[gene_sym]['syn']:
                print("\t\t{0}".format(syn.strip()))
                cursor.execute( add_gs, (gene_cache[gene_sym], syn.strip()) )

            print("\tLINKS:")
            for link in dstore[gene_sym]['links']:
                print("\t\t{0} - {1}".format(link['label'], link['href']))

                if link['label'].startswith('MGI'):
                    loaderutils.add_gene_url(cursor, gene_id, 'MGI', link['href'])
                elif link['label'].startswith('IGTC'):
                    loaderutils.add_gene_url(cursor, gene_id, 'IGTC', link['href'])
                elif link['label'].startswith('IMSR'):
                    loaderutils.add_gene_url(cursor, gene_id, 'IMSR', link['href'])
                elif link['label'].startswith('BioGPS'):
                    loaderutils.add_gene_url(cursor, gene_id, 'BioGPS', link['href'])
                elif link['label'] == 'GEO':
                    loaderutils.add_gene_url(cursor, gene_id, 'GEO', link['href'])
                elif 'phenome.jax.org' in link['href']:
                    continue
                    #loaderutils.add_gene_url(cursor, gene_id, 'MPD', link['href'])
    

    cnx.commit()
    cursor.close()
    cnx.close()
    
if __name__ == '__main__':
    main()







