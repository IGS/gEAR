#!/opt/bin/python3

"""
Saves a GeneCart object (from genecart.js) serialized as
formData and submitted to this script.  This assumes we're creating
a NEW GeneCart object
"""

import cgi, json
import re, sys
from pathlib import Path

TWO_LEVELS_UP = 2
abs_path_gear = Path(__file__).resolve().parents[TWO_LEVELS_UP]
abs_path_lib = abs_path_gear.joinpath('lib')
sys.path.insert(0, str(abs_path_lib))
import geardb
from gear.userhistory import UserHistory

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")


def validate_weighted_gene_cart(df):
    """Ensure weighted gene cart meets the requirements.  Returns a boolean."""
    # Check that first column is identifiers, second column is gene symbols, and following columns are numeric weights
    if len(df.columns) < 3:
        return False

    # Objects are variable string length.
    if df[df.columns[0]].dtype not in ["string", "object"] \
        or df[df.columns[1]].dtype not in ["string", "object"]:
        return False

    # The third column onward must be a numeric weight
    for col in df[df.columns[2:]]:
        if df[col].dtype != 'float64':
            return False

    # The first column has to be unique identifiers
    if not df[df.columns[0]].is_unique:
        return False

    return True

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
        pasted_genes = pasted_genes.replace('\n', ' ').replace('\r', '').replace('\t', ' ')

        for gene_sym in pasted_genes.split(' '):
            if len(gene_sym) > 0:
                gene = geardb.Gene(gene_symbol=gene_sym)
                gc.add_gene(gene)

    elif upload_type == 'uploaded-unweighted':
        gc.gctype = 'unweighted-list'

        fileitem = form['new_cart_file']
        if fileitem.filename:
            pasted_genes = fileitem.file.read().decode().replace(",", " ")
            pasted_genes = re.sub(r"\s+", " ", pasted_genes)
            pasted_genes = pasted_genes.replace('\n', ' ').replace('\r', '').replace('\t', ' ')

            for gene_sym in pasted_genes.split(' '):
                gene = geardb.Gene(gene_symbol=gene_sym)
                gc.add_gene(gene)
        else:
            raise Exception("Didn't detect an uploaded file for an uploaded-unweighted submission")

    elif upload_type == 'uploaded-weighted':
        import anndata
        import pandas as pd
        import string
        gc.gctype = 'weighted-list'

        # sanitize the file name
        fileitem = form['new_cart_file']
        valid_chars = "-_.() {}{}".format(string.ascii_letters, string.digits)
        fileitem.filename = ''.join(c for c in fileitem.filename if c in valid_chars)

        source_file_ext = ".tab"
        source_file_path = CARTS_BASE_DIR.joinpath("cart.{0}{1}".format(gc.share_id, source_file_ext))
        h5dest_file_path = CARTS_BASE_DIR.joinpath("cart.{0}.h5ad".format(gc.share_id))

        df = None
        try:
            if fileitem.filename.endswith('xlsx') or fileitem.filename.endswith('xls'):
                df = pd.read_excel(fileitem.file, sheet_name=0)
            elif fileitem.filename.endswith('tab') or fileitem.filename.endswith('tsv'):
                df = pd.read_csv(fileitem.file, sep='\t')
            elif fileitem.filename.endswith('csv'):
                df = pd.read_csv(fileitem.file, sep=',')
            else:
                raise Exception("Unsupported file type for carts uploaded. File name: {0}. Supported extensions: ['xlsx', 'xls', 'tab', 'tsv', 'csv']".format(fileitem.filename))

            is_valid = validate_weighted_gene_cart(df)

            if not is_valid:
                raise Exception("Weighted gene cart is not valid. Ensure first column is unique identifiers, second column is gene symbols, and following columns are numeric weights.")

            # Write dataframe to tab file
            try:
                df.to_csv(source_file_path, sep='\t', index=False)
            except:
                raise Exception("Could not write data to tab file: {0}".format(source_file_path))

            # First two columns make adata.var
            var = df[df.columns[:2]]
            var.set_index(var.columns[0], inplace=True)
            for gene_sym in var[var.columns[0]]:
                gene = geardb.Gene(gene_symbol=gene_sym)
                gc.add_gene(gene)

            # Remaining columns make adata.X
            X = df[df.columns[2:]].transpose().to_numpy()
            obs = pd.DataFrame(index=df.columns[2:])
            # Create the anndata object and write to h5ad
            adata = anndata.AnnData(X=X, obs=obs, var=var)
            adata.write(filename=h5dest_file_path)

        except Exception as e:
            print(str(e))
            sys.exit(1)
    elif upload_type == "labeled-list":
        raise NotImplementedError("Not implemented")
    else:
        raise Exception("Invalid upload type: {0}".format(upload_type))

    gc.save()

    result = { 'id': gc.id }
    print(json.dumps(result))

    if gc.user_id:
        history = UserHistory()
        history.add_record(
            user_id=gc.user_id,
            entry_category='gene_cart_added',
            label="Gene cart added: '{0}'".format(gc.label),
            gene_cart_share_id=gc.share_id
        )

if __name__ == '__main__':
    main()
