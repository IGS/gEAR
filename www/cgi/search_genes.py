#!/opt/bin/python3
"""
Searches the reference annotations for genes and returns them and their complete annotations.-
"""

import cgi, json
import sys
import mysql.connector
import os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
# TODO: next line will go away after conversion to geardb is complete
import gear.db

# results will not be shown after this count
MAX_GENE_SEARCH_LIMIT = 100000

def main():
    form = cgi.FieldStorage()

    ## can search for more than one gene symbol, separated by spaces
    search_gene_symbol = form.getvalue('search_gene_symbol')

    # Turn any commas into spaces
    search_gene_symbol = search_gene_symbol.replace(',', ' ')

    exact_match = form.getvalue('exact_match')

    # Get list of gene_ids found in miRNA_family table
    # TODO: refactor these to be on gene name rather than int (if they are)
    # TODO: check the previous revision to see how these were handled if miRNA
    #cached_mirna_family_gene_ids = get_mirna_family_gene_ids(cursor)

    gene_c = geardb.GeneCollection()
    gene_c.get_by_gene_symbol(gene_symbol=search_gene_symbol, exact=exact_match)

    result = { 'genes': gene_c.genes }

    for gene in result['genes']:
        gene.load_go_terms()
        gene.load_dbxref_links()
        gene.load_aliases()

        # these were in the previous result output.  unsure if they will be retained
        #gene['aliases_searched']
        #gene['mirna_family']

    syms = dict()
    for gene in result['genes']:
        ## has this symbol been seen?
        gene_sym = gene.gene_symbol.lower()

        if gene_sym not in syms:
            syms[gene_sym] = {'by_organism': dict()}

        if gene.organism_id not in syms[gene_sym]['by_organism']:
            syms[gene_sym]['by_organism'][gene.organism_id] = list()

        syms[gene_sym]['by_organism'][gene.organism_id].append(gene)

    print('Content-Type: application/json\n\n')
    print(json.dumps(syms))

def get_mirna_family_gene_ids(cursor):
    cached_mirna_ids = {}

    qry_stemloop = '''
                SELECT DISTINCT(stem_loop_id)
                FROM mirna_family
                '''
    cursor.execute(qry_stemloop,)
    for row in cursor:
        cached_mirna_ids[ row[0] ] = 'stem-loop'

    qry_mature = '''
                SELECT DISTINCT(mature_id)
                FROM mirna_family
                '''
    cursor.execute(qry_mature,)
    for row in cursor:
        cached_mirna_ids[ row[0] ] = 'mature'

    return cached_mirna_ids

if __name__ == '__main__':
    main()
