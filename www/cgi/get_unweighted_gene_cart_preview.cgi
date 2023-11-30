#!/opt/bin/python3

"""

"""

import cgi
import json
import os
import string
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    share_id = form.getvalue('share_id')
    valid_chars = "%s%s" % (string.ascii_letters, string.digits)
    share_id = ''.join(c for c in share_id if c in valid_chars)

    result = { 'success': 0 }

    # Get genes from the gene cart
    if not share_id:
        result['message'] = 'No share ID provided'
        result["success"] = 0
        print(json.dumps(result))
        return

    gene_cart = geardb.get_gene_cart_by_share_id(share_id)
    if not gene_cart:
        result['message'] = 'Invalid share ID'
        result["success"] = 0
        print(json.dumps(result))
        return

    genes = gene_cart.genes
    organism_id = gene_cart.organism_id

    # Get the Gene objects for this gene cart
    gene_collection = geardb.GeneCollection()
    gene_collection.get_by_gene_symbol(gene_symbol=' '.join(genes), organism_id=organism_id, exact=True)


    gene_dict = {} # Key is ensembl ID, value is dict of gene symbol and product.
    # Create an dict of gene symbol to gene object
    for gene in gene_collection.genes:
        gene_dict[gene.ensembl_id] = { 'gene_symbol': gene.gene_symbol, 'product': gene.product }

    result["gene_info"] = gene_dict
    result['success'] = 1
    print(json.dumps(result))

if __name__ == '__main__':
    main()
