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

def get_standard_images(cursor=None, gene_id=None):
    ## TODO: Make this dynamic later:
    image_standard_datasets = ['liasdf97-e9a2-po1u-kj11-1k282bjg8j81', 'liasdf96-oi12-n132-i8j2-9sd82bjg8iu9', 'liasdf95-b723-n132-02kd-2dd8yujg8jsi', '9ij23oiu-l12n-oisn-12bn-123b8b0982bn']
    images = list()

    for dataset_id in image_standard_datasets:
        if os.path.isfile("../img/standard/{0}/{1}.png".format(dataset_id, gene_id)):
            images.append({'label': dataset_id,
                            'ensembl_id': gene_id,
                            'image_url': "./img/standard/{0}/{1}.png".format(dataset_id, gene_id)})
        else:
            images.append({'label': dataset_id,
                            'ensembl_id': gene_id,
                            'image_url': "".format(dataset_id, gene_id)})

    return images


def get_supplemental_images(cursor=None, label=None, symbol=None):
    qry = "SELECT ensembl_id, image_url FROM supplemental_images WHERE label = %s AND gene_symbol = %s LIMIT 1"
    cursor.execute(qry, (label, symbol))

    for row in cursor:
        return [{'label': label, 'ensembl_id': row[0], 'image_url': row[1]}]

if __name__ == '__main__':
    main()
