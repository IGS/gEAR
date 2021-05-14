#!/opt/bin/python3

"""

"""

import cgi, json
import os, sys

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

# this is needed so that we don't get TclError failures in the underlying modules
import matplotlib
matplotlib.use('Agg')

import numpy as np

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

    genes_prefix = form.getvalue('genes_prefix')
    filter_mito_perc = form.getvalue('filter_mito_perc')
    filter_mito_count = form.getvalue('filter_mito_count')
    save_dataset = int(form.getvalue('save_dataset'))

    adata = ana.get_adata()

    # primary or public analyses won't be after this
    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path()
    dest_directory = os.path.dirname(dest_datafile_path)

    if not os.path.exists(dest_directory):
        os.makedirs(dest_directory)

    # Case-insensitive check.  series.str.function returns a series, so add str again
    mito_genes = adata.var.gene_symbol.str.lower().str.startswith(genes_prefix.lower())

    # for each cell compute fraction of counts in mito genes vs. all genes
    # the ".A1" is only necessary, as X is sparse - it transform to a dense array after summing
    try:
        adata.obs['percent_mito'] = np.sum(
            adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

        # add the total counts per cell as observations-annotation to adata
        adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1

    except AttributeError:
        adata.obs['percent_mito'] = np.sum(
            adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)

        # add the total counts per cell as observations-annotation to adata
        adata.obs['n_counts'] = np.sum(adata.X, axis=1)

    os.chdir(os.path.dirname(dest_datafile_path))
    axs = sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
                       jitter=0.4, multi_panel=True, save="_qc_by_mito.png")

    ax = sc.pl.scatter(adata, x='n_counts', y='percent_mito', save="_percent_mito.png")
    ax = sc.pl.scatter(adata, x='n_counts', y='n_genes', save="_n_genes.png")

    result = {'success': 1}

    if save_dataset:
        if filter_mito_count:
            adata = adata[adata.obs['n_genes'] < int(filter_mito_count), :]

        if filter_mito_perc:
            adata = adata[adata.obs['percent_mito'] < float(filter_mito_perc), :]

        adata.write(dest_datafile_path)
        (n_obs, n_genes) = adata.shape
        result['n_obs'] = n_obs
        result['n_genes'] = n_genes

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
