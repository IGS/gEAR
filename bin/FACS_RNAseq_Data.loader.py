#!/usr/bin/env python3

"""

First run this, then don't forget to run the annotation loader from the GBK files so that
product names and GO terms get added.

BEFORE:
    mysql> select is_primary, count(*) from gene_symbol group by is_primary;
    +------------+----------+
    | is_primary | count(*) |
    +------------+----------+
    |          0 |     2602 |
    |          1 |   142200 |
    +------------+----------+

    mysql> select count(*) from gene;
    +----------+
    | count(*) |
    +----------+
    |   142200 |
    +----------+

AFTER:
    mysql> select is_primary, count(*) from gene_symbol group by is_primary;
    +------------+----------+
    | is_primary | count(*) |
    +------------+----------+
    |          0 |    40123 |
    |          1 |   159156 |
    +------------+----------+

    mysql> select count(*) from gene;
    +----------+
    | count(*) |
    +----------+
    |   159156 |
    +----------+

File: FACS_RNAseq_Data.tab
 Columns:
 0 - MGI ID (ex: MGI:101757)
 1 - Symbol
 2 - mm9_coordinates
 3 - Chr	
 4 - cM Position	
 5 - mm10_start	
 6 - mm10_end	
 7 - strand	
 8 - Description	
 9 - Type	
10 - Synonyms (pipe separated)
11 - old gene name
....

File: MGI_Gene_Model_Coord.rpt
 Columns: 
 0. MGI accession id     
 1. marker type  
 2. marker symbol        
 3. marker name  
 4. genome build 
 5. Entrez gene id       
 6. NCBI gene chromosome 
 7. NCBI gene start      
 8. NCBI gene end        
 9. NCBI gene strand    
10. Ensembl gene id     
11. Ensembl gene chromosome     
12. Ensembl gene start  
13. Ensembl gene end    
14. Ensembl gene strand 
15. VEGA gene id        
16. VEGA gene chromosome        
17. VEGA gene start     
18. VEGA gene end       
19. VEGA gene strand

"""

import mysql.connector
import configparser
import loaderutils
import os
import re

def main():
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
    gene_ensembl_ids = loaderutils.cache_genes_by_ensembl_id(cursor)
    gene_aliases = loaderutils.cache_gene_aliases(cursor, True)
    mgi_aliases = dict()

    gear_bin_path = os.path.dirname(os.path.realpath(__file__))
    gear_root_path = os.path.dirname(gear_bin_path)
    
    for line in open('{}/datasets/FACS_RNAseq_Data.tab'.format(gear_root_path)):
        line = line.rstrip()
        if line.startswith('MGI:'):
            cols = line.split("\t")
            mgi_id = cols[0]
            gene_symbol = cols[1]
            mgi_aliases[mgi_id] = list()
            mgi_aliases[mgi_id].append(gene_symbol)

            aliases = cols[10].split('|')
            for alias in aliases:
                if ' ' not in alias and len(alias) > 0:
                    mgi_aliases[mgi_id].append(alias)

            #input("Symbol: {0} has {1} aliases: {2}".format(gene_symbol, len(mgi_aliases[mgi_id]), mgi_aliases[mgi_id]))

    for line in open('{}/datasets/MGI_Gene_Model_Coord.rpt'.format(gear_root_path)):
        line = line.rstrip()
        if line.startswith('MGI'):
            cols = line.split("\t")
            mgi_id = cols[0]
            ensembl_id = cols[10]
            gene_symbol = cols[2]
            lc_gene_symbol = cols[2].lower()

            if ensembl_id == 'null':
                continue
            
            print("Checking for ENSEMBL_ID: ({0}) ... ".format(ensembl_id), end='')

            if ensembl_id in gene_ensembl_ids:
                gene_id = gene_ensembl_ids[ensembl_id]
                # This gene is already present.  Add any secondary gene symbols
                print(" found")
                if lc_gene_symbol != gene_aliases[gene_id]['p']:
                    print("Primary gene symbol ({0}) for gene_id:{1} didn't match MGI file value of {2}".format(gene_aliases[gene_id]['p'], gene_id, gene_symbol))
                    # Add it as a secondary gene_symbol, if it isn't already one
                    if lc_gene_symbol not in gene_aliases[gene_id]['s'] and len(lc_gene_symbol) <= 20:
                        loaderutils.add_gene_symbol(cursor, gene_id, gene_symbol, 0)
                        gene_aliases[gene_id]['s'].append(gene_symbol.lower())
            else:
                # Add the gene.  This automatically does a primary gene symbol too
                print(" new!")
                gene_id = loaderutils.add_gene(cursor, ensembl_id, 1, gene_symbol, 'protein_coding')
                gene_aliases[gene_id] = {'p': gene_symbol.lower(), 's': []}

            # process any SHIELD aliases
            if mgi_id in mgi_aliases:
                for alias in mgi_aliases[mgi_id]:
                    lc_alias = alias.lower()
                    if lc_alias != gene_aliases[gene_id]['p'] and lc_alias not in gene_aliases[gene_id]['s'] and len(lc_alias) <= 20:
                        loaderutils.add_gene_symbol(cursor, gene_id, alias, 0)
                        gene_aliases[gene_id]['s'].append(lc_alias)

    cursor.close()
    cnx.commit()
    cnx.close()

if __name__ == '__main__':
    main()


