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
    
    ana = geardb.Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id, session_id=session_id )

    filter_cells_lt_n_genes = form.getvalue('filter_cells_lt_n_genes')
    filter_cells_gt_n_genes = form.getvalue('filter_cells_gt_n_genes')
    filter_genes_lt_n_cells = form.getvalue('filter_genes_lt_n_cells')
    filter_genes_gt_n_cells = form.getvalue('filter_genes_gt_n_cells')

    adata = ana.get_adata()
    dest_datafile_path = None

    if ana.type == 'primary':
        ana.type = 'user_unsaved'
        dest_datafile_path = ana.dataset_path()
    else:
        raise Exception("ERROR: unsupported analysis type: {0}".format(ana.type))

    dest_directory = os.path.dirname(dest_datafile_path)

    if not os.path.exists(dest_directory):
        os.makedirs(dest_directory)

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # API documentation states one filter param per call
    if filter_genes_lt_n_cells:
        sc.pp.filter_genes(adata, min_cells=int(filter_genes_lt_n_cells))

    if filter_genes_gt_n_cells:
        sc.pp.filter_genes(adata, max_cells=int(filter_genes_gt_n_cells))
        
    if filter_cells_lt_n_genes:
        sc.pp.filter_cells(adata, min_genes=int(filter_cells_lt_n_genes))

    if filter_cells_gt_n_genes:
        sc.pp.filter_cells(adata, max_genes=int(filter_cells_gt_n_genes))

    adata.write(dest_datafile_path)
    (n_obs, n_genes) = adata.shape
    ## I don't see how to get the save options to specify a directory
    # sc.settings.figdir = 'whateverpathyoulike' # scanpy issue #73
    os.chdir(os.path.dirname(dest_datafile_path))

    ax = sc.pl.highest_expr_genes(adata, n_top=20, gene_symbols='gene_symbol', save=".png")

    result = {'n_obs': n_obs, 'n_genes': n_genes}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
