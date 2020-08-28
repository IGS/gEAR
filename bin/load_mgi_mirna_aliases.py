#!/opt/bin/python3

'''
MGI entries are tab-delimited. Like this:

MGI:2676793	Mirlet7a-1	O	microRNA let7a-1	24.84	13	Gene		387244	Mirnlet7a|mmu-let-7a-1|Mirnlet7a-1|Let-7a	miRNA gene	48538179	48538272	-	miRNA|ncRNA


'''

import argparse
import configparser
import cgi, json
import mysql.connector
from mysql.connector import errorcode
import re
import os
import sys

import csv

sys.path.append("{0}/../".format(os.path.dirname(sys.argv[0])))
import loaderutils

def main():
    parser = argparse.ArgumentParser( description='Process MGI miR and miRBase symbols')
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to MGI_EntrezGene file' )
    args = parser.parse_args()

    # This will hold all miRNA info
    miRNAs = {}

    print('\nConnecting to gEAR database...\n')
    config = configparser.ConfigParser()
    config.read('../gear.ini')

    cnx = None
    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'])
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password", file=sys.stderr)
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist", file=sys.stderr)
        else:
            print(err, file=sys.stderr)

    print("  ...Done.\n")

    cursor = cnx.cursor()

    print("Caching genes from gEAR...")
    # Get aliases from gEAR database
    cached_genes = loaderutils.cache_genes_by_gene_symbol(cursor, lower=False)
    cached_primary_genes = loaderutils.cache_genes_by_primary_gene_symbol(cursor, lower=False)
    print("gEAR genes have been cached.\n")

    print("Parsing MGI_EntrezGene file and adding aliases to gEAR portal...\n")
    mgi_file = open(args.file)
    reader = csv.reader(mgi_file, delimiter='\t')

    counter_added = 0
    counter_changed = 0
    for cols in reader:
        # print(cols)
        if cols[0].startswith('MGI:'):
            mgi_gene_sym = cols[1]
            # mgi_gene_sym = cols[1].lower() #column: Marker Symbol
            mgi_gene_name = cols[3] #column: Marker Name

            #Only interested in miR entries
            if mgi_gene_sym.startswith( ('mir', 'Mir') ) and mgi_gene_name.startswith('microRNA'):
                    # print(cols)
                    # Skip MGI miR not found in gEAR
                    if mgi_gene_sym not in cached_genes and mgi_gene_sym.lower() not in cached_genes:
                        # if mgi_gene_sym == 'Mir6415':
                        #     print("{0} is being skipped. {1}".format(mgi_gene_sym, cols[9]) )
                        continue
                    else:
                        mgi_aliases = cols[9] #column: 'Synonyms |-delimited'
                        mgi_aliases_list = []
                        # Multiple aliases are '|' separated
                        if '|' in mgi_aliases:
                            mgi_aliases_list = mgi_aliases.split('|')
                        else:
                            mgi_aliases_list.append(mgi_aliases)

                        # TODO: why is this miR skipped in processing?!
                        # if 'mmu-mir-6415' in mgi_aliases_list:
                        #     print('alias: {0} LINE 85'.format(alias))

                        if mgi_gene_sym in cached_primary_genes or mgi_gene_sym.lower() in cached_primary_genes:
                            # MGI sym is primary. Add any aliases to gene_symbol table
                            for alias in mgi_aliases_list:
                                # Only interested in mouse miRNAs.
                                if alias.startswith('mmu'):
                                    # if alias == 'mmu-mir-6415':
                                    #     print('alias: {0} LINE 93'.format(alias))
                                    result = process_alias(cnx, cursor, mgi_gene_sym, alias, cached_genes, cached_primary_genes)

                                    if result == 'alias_changed':
                                        counter_changed += 1

                                    if result == 'alias_added':
                                        counter_added += 1

                        else:
                            # MGI sym is not primary. Change it to primary
                            try:
                                mgi_gene_id = cached_genes[mgi_gene_sym]
                            except KeyError:
                                mgi_gene_id = cached_genes[mgi_gene_sym.lower()]

                            make_sym_primary = "UPDATE gene_symbol SET is_primary = %s WHERE gene_id = %s AND label = %s"
                            cursor.execute( make_sym_primary, (1, mgi_gene_id, mgi_gene_sym) )
                            cnx.commit()
                            print("{0} is now a primary symbol".format(mgi_gene_sym))

                            # If an alias listed is marked 'is_primary = 1' change it to '0'
                            for alias in mgi_aliases_list:
                                # Only interested in mouse miRNAs.
                                if alias.startswith('mmu'):
                                    # if alias == 'mmu-mir-6415':
                                    #     print('alias: {0} LINE: 119'.format(alias))
                                    result = process_alias(cnx, cursor, mgi_gene_sym, alias, cached_genes, cached_primary_genes)

                                    if result == 'alias_changed':
                                        counter_changed += 1

                                    if result == 'alias_added':
                                        counter_added += 1

    cursor.close()
    cnx.close()
    print("\n\nFinished\n")
    print("Total number of miR aliases added to gEAR: {0}".format(counter_added))
    print("Total number of miR aliases were mistakenly marked as primary symbols: {0}".format(counter_changed))


def process_alias(cnx, cursor, mgi_gene_sym, alias, cached_genes, cached_primary_genes):
    # Is alias in gEAR?
    if alias in cached_genes:
        # Is alias marked as primary sym? (by mistake)
        if alias in cached_primary_genes:
            alias_id = cached_primary_genes[alias]

            # alias is mistakenly set as primary - change it to is_primary=0
            make_sym_secondary = "UPDATE gene_symbol SET is_primary = %s WHERE gene_id = %s AND label = %s"
            cursor.execute( make_sym_secondary, (0, alias_id, alias) )
            cnx.commit()
            print("{0} changed to secondary".format(alias))

            return 'alias_changed'

    #alias NOT in gEAR
    else:
        # Add alias (marked as secondary sym)
        try:
            mgi_gene_id = cached_genes[mgi_gene_sym]
        except KeyError:
            mgi_gene_id = cached_genes[mgi_gene_sym.lower()]
        loaderutils.add_gene_symbol(cursor, mgi_gene_id, alias, 0)
        cnx.commit()
        print("{0} was added as secondary symbol".format(alias))

        return 'alias_added'


if __name__ == '__main__':
    main()
