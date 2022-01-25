#!/opt/bin/python3

"""
Saves a GeneCart object (from genecart.js) serialized as
formData and submitted to this script.  This assumes we're creating
a NEW GeneCart object
"""

import cgi, json
import os, re, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
from pprint import pprint
import geardb

def main():
    print('Content-Type: application/json\n\n')
    gc = geardb.GeneCart()
    form = cgi.FieldStorage()

    gc.label = form.getvalue('new_cart_label')
    gc.organism_id = form.getvalue('new_cart_organism_id')
    gc.ldesc = form.getvalue('new_cart_ldesc')
    gc.is_public = form.getvalue('is_public')

    user_logged_in = geardb.get_user_from_session_id(form.getvalue('session_id'))
    gc.user_id = user_logged_in.id

    upload_type = form.getvalue('new_cart_upload_type')

    if upload_type == 'pasted_genes':
        gc.gctype = 'unweighted-list'
        pasted_genes = form.getvalue('new_cart_pasted_genes').replace(',', ' ').replace('  ', ' ')

        for gene_sym in pasted_genes.split(' '):
            gene = geardb.Gene(gene_symbol=gene_sym)
            gc.add_gene(gene)

    elif upload_type == 'uploaded-unweighted':
        gc.gctype = 'unweighted-list'

        fileitem = form['new_cart_file']
        if fileitem.filename:
            pasted_genes = fileitem.file.read().decode().replace(",", " ")
            pasted_genes = re.sub(r"\s+", " ", pasted_genes)

            for gene_sym in pasted_genes.split(' '):
                gene = geardb.Gene(gene_symbol=gene_sym)
                gc.add_gene(gene)
        else:
            raise Exception("Didn't detect an uploaded file for an uploaded-unweighted submission")

    elif upload_type == 'uploaded-weighted':
        import scanpy as sc
        import string
        gc.gctype = 'weighted-list'

        # sanitize the file name
        fileitem = form['new_cart_file']
        valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
        fileitem.filename = ''.join(c for c in fileitem.filename if c in valid_chars)

        file_ext = os.path.splitext(fileitem.filename)[1]
        package_dir = os.path.dirname(os.path.abspath(__file__))
        carts_dir =  os.path.join(package_dir, '..', 'carts')
        source_file_path = os.path.join(carts_dir, "cart.{0}{1}".format(gc.share_id, file_ext))
        h5dest_file_path = os.path.join(carts_dir, "cart.{0}.h5ad".format(gc.share_id))

        with open(source_file_path, 'wb') as sfh:
            sfh.write(fileitem.file.read())

        if fileitem.filename.endswith('xlsx') or fileitem.filename.endswith('xls'):
            adata = sc.read_excel(source_file_path, index_col=0).transpose()
            adata.write(filename=h5dest_file_path)

        elif fileitem.filename.endswith('tab'):
            adata = sc.read_csv(source_file_path, delimiter="\t", first_column_names=True).transpose()
            adata.write(filename=h5dest_file_path)

        elif fileitem.filename.endswith('csv'):
            adata = sc.read_csv(source_file_path, first_column_names=True).transpose()
            adata.write(filename=h5dest_file_path)

        else:
            raise Exception("Unsupported file type for carts uploaded. File name: {0}".format(fileitem.filename))

    gc.save()

    result = { 'id': gc.id }
    print(json.dumps(result))

if __name__ == '__main__':
    main()
