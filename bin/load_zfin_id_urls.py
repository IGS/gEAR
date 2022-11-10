#!/opt/bin/python3

'''
Scan ZFin id list and parse out ZDB-GENE ids and coorresponding external source ids

2 Examples:
ZDB-GENE-000112-32	U93460,ZDB-GENE-000112-32,couptf3
ZDB-GENE-030909-1	378723,5379,AAH56743,AAH65656,AB242329,BAE48583,BC056743,BC065656,ENSDARG00000070913,ENSDARP00000095266,IPR009071,IPR022097,IPR032643,NM_213118,NP_998283,PF00505,PF12336,PS50118,Q6P0E1,ZDB-GENE-030131-1891,ZDB-GENE-030131-2647,ZDB-GENE-030909-1,ZDB-GENE-040426-2124,ZDB-GENE-040426-2425,sox2

Notes:

Skip all except: lines starting with "ZDB-GENE-"

'''

import argparse
import mysql.connector
import configparser
import os
import re
import loaderutils

def main():
    parser = argparse.ArgumentParser( description='')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read('gear.ini')

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

    # query db to get list of gene_symbols & ids
    gene_cache = loaderutils.cache_genes_by_gene_symbol(cursor, lower=True)

    #print(gene_cache)

    current = dict()

    count = 0
    for line in open(args.input_file):
        line = line.rstrip()

        if line.startswith('ZDB-GENE-'):

            if 'gene_id' in current:
                # exclude entries that lack a gene_id
                if current['gene_id'] is not None:
                    print("Adding {0} - {1}".format(current['gene_id'], current['url']))

                    #add to db
                    add_zfin_url(cursor, current)
                    count += 1

                    #reset for next term
                    current = {'gene_id': None, 'label': 'ZFIN', 'url': None}

            zfin_pair = line.split()

            zfin_id = zfin_pair[0]
            symbols = zfin_pair[1].split(',') #split into list of symbols associated with ZFIN id

            for symbol in symbols:
                symbol = symbol.lower()
                if symbol in gene_cache:
                    current['gene_id'] = gene_cache[symbol]
                    current['label'] = 'ZFIN'
                    current['url'] = 'http://zfin.org/' + zfin_id

    print("Added {0} new urls.".format(count))

    cnx.commit()
    cursor.close()
    cnx.close()

def add_zfin_url( curs, atts ):
    add_term_sql = ("INSERT INTO gene_urls (gene_id, label, url) "
                    "VALUES (%s, %s, %s)")


    curs.execute(add_term_sql, (atts['gene_id'], atts['label'], atts['url']))

if __name__ == '__main__':
    main()







