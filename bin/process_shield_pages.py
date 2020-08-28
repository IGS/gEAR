#!/usr/bin/env python3

"""

"""

import mysql.connector
import configparser
import loaderutils
import os
import re

def main():
    INPUT_BASE_DIR = '/usr/local/projects/gEAR/shield_pages'
    URL_BASE = 'https://shield.hms.harvard.edu/viewgene.html?gene='

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

    gene_syms = loaderutils.cache_genes_by_primary_gene_symbol(cursor, lower=True)

    for file in os.listdir(INPUT_BASE_DIR):
        file_path = "{0}/{1}".format(INPUT_BASE_DIR, file)
        
        m = re.match(".+\?gene=(.+)", file)
        if m:
            gene_id = m.group(1)
        else:
            continue

        lower_gene_id = gene_id.lower()

        if lower_gene_id not in gene_syms:
            print("PATH: processing file for UNKNOWN gene: {0}".format(gene_id))

        aliases = list()
        parsing_synonyms = False
        image_url = None

        for line in open(file_path):
            if parsing_synonyms == True:
                m = re.match("(.+)<\/td>", line)
                if m:
                    aliases = m.group(1).split("; ")
                elif '</tr>' in line:
                    parsing_synonyms = False

                continue
                
            if 'id="gene_synonyms"' in line:
                parsing_synonyms = True
                continue
            elif 'id="FACS_chart"' in line:
                m = re.search("href='(http.+.png)'", line)
                if m:
                    image_url = m.group(1)
                else:
                    raise Exception("ERROR: failed to find a URL for gene {0}".format(gene_id))

        novel_aliases = list()
        #for alias in aliases:
        #    if alias.lower() not in gene_syms[lower_gene_id]:
        #        novel_aliases.append(alias)

        # ALTER TABLE gene add shield_facs_chart_url VARCHAR(255);
        print("gene_id:{4}\tgene_sym:{0}\turl:{1}\taliases:{3}/{2}".format(gene_id, image_url, len(novel_aliases), len(aliases), gene_syms[lower_gene_id]))
        if image_url is None:
            print("WARN: No FACS chart URL for gene: {0}".format(lower_gene_id))
            print("UPDATE gene SET shield_facs_chart_url = '' WHERE id = {0};".format(gene_syms[lower_gene_id]))
        else:
            print("UPDATE gene SET shield_facs_chart_url = '{0}' WHERE id = {1};".format(image_url, gene_syms[lower_gene_id]))
            
    cursor.close()
    cnx.close()


def run_command(cmd):
    #print("INFO: [{1}] Running command: {0}\n".format(cmd, datetime.datetime.now()), flush=True)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
       raise Exception("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))

if __name__ == '__main__':
    main()


