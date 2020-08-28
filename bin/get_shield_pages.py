#!/usr/bin/env python3

"""

Queries all the primary gene symbols from our internal database and grabs the corresponding
page from SHIELD, if present.  Example:

https://shield.hms.harvard.edu/viewgene.html?gene=Rfx3

Parsing and actually doing anything with these is the role of other scripts.

"""

import argparse
import mysql.connector
import configparser
import os
import re
import subprocess 

def main():
    parser = argparse.ArgumentParser( description='')

    OUTPUT_BASE_DIR = '/usr/local/projects/gEAR/shield_pages'
    URL_BASE = 'https://shield.hms.harvard.edu/viewgene.html?gene='

    ## output file to be written
    #parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
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

    os.chdir(OUTPUT_BASE_DIR)

    qry = ("SELECT DISTINCT label FROM gene_symbol WHERE is_primary = 1")
    search_terms = list()
    
    cursor.execute(qry, search_terms)

    for (label,) in cursor:
        # don't do it if it contains a space
        if ' ' in label:
            continue

        if '(' in label:
            continue

        # skip it if the file exists already
        if os.path.exists("viewgene.html?gene={0}".format(label)):
            continue
        
        
        print("Processing: {0}".format(label))
        url = "{0}{1}".format(URL_BASE, label)
        run_command("wget {0}".format(url))

    cursor.close()
    cnx.close()


def run_command(cmd):
    #print("INFO: [{1}] Running command: {0}\n".format(cmd, datetime.datetime.now()), flush=True)
    return subprocess.check_output(cmd, shell=True)
    



if __name__ == '__main__':
    main()


