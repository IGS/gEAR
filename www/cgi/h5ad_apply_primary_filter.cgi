#!/opt/bin/python3

"""
Used to apply filters to an H5AD dataset.  These include:

- Cells with < N genes
- Cells with > N genes
- Genes in < N cells
- Genes in > N cells

If the input dataset is a primary one, this script makes a copy
into a user/session specific directory first.

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

import scanpy as sc
sc.settings.verbosity = 0

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    result = {'success': 0, 'n_obs': None, 'n_genes': None}

    ana = geardb.Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id, session_id=session_id )

    filter_cells_lt_n_genes = form.getvalue('filter_cells_lt_n_genes')
    filter_cells_gt_n_genes = form.getvalue('filter_cells_gt_n_genes')
    filter_genes_lt_n_cells = form.getvalue('filter_genes_lt_n_cells')
    filter_genes_gt_n_cells = form.getvalue('filter_genes_gt_n_cells')

    adata = ana.get_adata()

    # This step should only be performed on the original dataset.
    # However the filtering will be saved in a temp directory to avoid overwriting the original dataset.
    if ana.type == 'primary':
        ana.type = 'user_unsaved'
    else:
        raise Exception("ERROR: unsupported analysis type: {0}".format(ana.type))

    dest_datafile_path = ana.dataset_path()
    dest_directory = os.path.dirname(dest_datafile_path)

    if not os.path.exists(dest_directory):
        os.makedirs(dest_directory)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    was_filtered = False

    # API documentation states one filter param per call
    if filter_genes_lt_n_cells:
        was_filtered = True
        sc.pp.filter_genes(adata, min_cells=int(filter_genes_lt_n_cells))

    if filter_genes_gt_n_cells:
        was_filtered = True
        sc.pp.filter_genes(adata, max_cells=int(filter_genes_gt_n_cells))

    if filter_cells_lt_n_genes:
        was_filtered = True
        sc.pp.filter_cells(adata, min_genes=int(filter_cells_lt_n_genes))

    if filter_cells_gt_n_genes:
        was_filtered = True
        sc.pp.filter_cells(adata, max_genes=int(filter_cells_gt_n_genes))

    # If no filters were selected, use initial dataset.
    # Filter to get the n_cells and n_genes obs columns
    if not was_filtered:
        sc.pp.filter_cells(adata, min_genes=0)

    adata.write(dest_datafile_path)
    (n_obs, n_genes) = adata.shape

    sc.settings.figdir = dest_directory + "/figures"

    try:
        ax = sc.pl.highest_expr_genes(adata, n_top=20, gene_symbols='gene_symbol', save=".png")
        result['success'] = 1
    except:
        result['success'] = 0

    result['n_obs'] = n_obs
    result['n_genes'] = n_genes

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
