#!/opt/bin/python3

"""

"""

import cgi, json
import os, sys, re

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

# this is needed so that we don't get TclError failures in the underlying modules
import matplotlib
matplotlib.use('Agg')

import scanpy as sc
sc.settings.verbosity = 0

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)

    ana = geardb.Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id,
                          session_id=session_id, user_id=user.id)

    genes_to_color = form.getvalue('genes_to_color')
    compute_pca = form.getvalue('compute_pca')

    if genes_to_color:
        genes_to_color = genes_to_color.replace(' ', '')

    if genes_to_color and ',' in genes_to_color:
        genes_to_color = genes_to_color.split(',')
    else:
        genes_to_color = [genes_to_color]

    adata = ana.get_adata()

    # primary or public analyses won't be after this
    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path()
    dest_directory = os.path.dirname(dest_datafile_path)

    if not os.path.exists(dest_directory):
        os.makedirs(dest_directory)

    # Compute PCA and make a scatter plot.
    if compute_pca == 'true':
        sc.tl.pca(adata, svd_solver='arpack')
        adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat R
        adata.write(dest_datafile_path)

    ## I don't see how to get the save options to specify a directory
    # sc.settings.figdir = 'whateverpathyoulike' # scanpy issue #73
    os.chdir(os.path.dirname(dest_datafile_path))

    missing_gene = None

    if genes_to_color:
        # Catch the error if any gene names are passed which aren't in the dataset
        try:
            # This can error like: ValueError: "RFX7" is invalid! specify valid sample annotation
            # https://scanpy.readthedocs.io/en/latest/api/scanpy.pl.pca.html#scanpy.pl.pca
            adata.var = adata.var.reset_index().set_index('gene_symbol')

            # ScanPy's pca function checks if there's a raw attribute
            # within adata. If raw attribute exists and use_raw argument is None,
            # pca decides to set use_raw to true and uses that. So, we need to
            # reinitialize Raw with our new adata. So it searches for gene symbols
            # in the newer index.
            # https://github.com/theislab/anndata/blob/master/anndata/base.py#L1020-L1022
            if adata.raw is not None:
                adata.raw = adata

            ax = sc.pl.pca(adata, color=genes_to_color, save=".png")
        except ValueError as err:
            m = re.search("\: (.+?) is not a valid", str(err))
            if m:
                missing_gene = m.groups(1)
            else:
                missing_gene = 'Unknown'
    else:
        ax = sc.pl.pca(adata, right_margin=0.2, save=".png")

    if missing_gene is None:
        sc.pl.pca_variance_ratio(adata, log=True, save=".png")
        result = {'success': 1, 'missing_gene': ''}
    else:
        result = {'success': 0, 'missing_gene': missing_gene}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

