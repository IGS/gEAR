#!/opt/bin/python3
"""
Searches the reference annotations for genes and returns them and their complete annotations.-
"""

import cgi, json
import sys
import os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

from gear.userhistory import UserHistory

# results will not be shown after this count
MAX_GENE_SEARCH_LIMIT = 100000

def main():
    print("DEBUG: search_genes.cgi logged", file=sys.stderr)
    form = cgi.FieldStorage()

    ## can search for more than one gene symbol, separated by spaces
    search_gene_symbol = form.getvalue('search_gene_symbol')

    # Turn any commas into spaces
    search_gene_symbol = search_gene_symbol.replace(',', ' ')

    exact_match = form.getvalue('exact_match')
    layout_share_id = form.getvalue('layout_share_id')
    user_id = form.getvalue('user_id')

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

    # log the search if user info is available
    if user_id:
        layout = geardb.get_layout_by_share_id(layout_share_id)

        ## shorten the gene string if it's silly long
        if len(search_gene_symbol) > 10:
            gene_symbol_label = "{0} ...".format(search_gene_symbol[0:50])
        else:
            gene_symbol_label = search_gene_symbol

        print("DEBUG: Adding history record", file=sys.stderr)
        history = UserHistory()
        history.add_record(
            user_id=user_id,
            entry_category='gene_search',
            label="\"{0}\" in {1}".format(gene_symbol_label, layout.label),
            gene_symbol=search_gene_symbol,
            layout_share_id=layout_share_id
        )

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
