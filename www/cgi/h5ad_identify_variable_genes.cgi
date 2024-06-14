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

    norm_counts_per_cell = float(form.getvalue('norm_counts_per_cell'))
    flavor = form.getvalue('flavor')
    n_top_genes = form.getvalue('n_top_genes', None)
    min_mean = float(form.getvalue('min_mean'))
    max_mean = float(form.getvalue('max_mean'))
    min_dispersion = float(form.getvalue('min_dispersion'))
    regress_out = form.getvalue('regress_out')
    scale_unit_variance = form.getvalue('scale_unit_variance')
    save_dataset = int(form.getvalue('save_dataset'))

    adata = ana.get_adata()

    if n_top_genes:
        n_top_genes = int(n_top_genes)

    # primary or public analysis should not be overwritten
    # this will alter the analysis object save destination
    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path()
    dest_directory = os.path.dirname(dest_datafile_path)

    if not os.path.exists(dest_directory):
        os.makedirs(dest_directory)

    # Per-cell normalize the data matrix, identify highly-variable genes and compute logarithm.
    sc.pp.normalize_total(adata, target_sum=norm_counts_per_cell)

    # log the data
    sc.pp.log1p(adata)
    adata.raw = adata

    # SAdkins - For some reason, if adata is backed these will fail when a copy of adata.X is made by the function.
    if n_top_genes:
        sc.pp.highly_variable_genes(
            adata, flavor=flavor, n_top_genes=n_top_genes)
    else:
        sc.pp.highly_variable_genes(
            adata, flavor=flavor, min_mean=min_mean, max_mean=max_mean, min_disp=min_dispersion)

    ## I don't see how to get the save options to specify a directory
    os.chdir(os.path.dirname(dest_datafile_path))
    sc.pl.highly_variable_genes(adata, save=".png")

    if save_dataset:
        # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed.
        if regress_out == 'true':
            sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

        if scale_unit_variance == 'true':
            sc.pp.scale(adata, max_value=10)

        adata.write(dest_datafile_path)

    ## actually do the filter so we can get the shape
    adata = adata[:, adata.var['highly_variable']]

    (n_obs, n_genes) = adata.shape

    # get a list of the top highly variable genes to display, sorted by best normalized dispersion
    highly_variable_genes = adata.var[adata.var.highly_variable].sort_values('dispersions_norm', ascending=False).gene_symbol
    if len(highly_variable_genes) > 5:
        highly_variable_genes = highly_variable_genes[:5]

    top_genes = ", ".join(highly_variable_genes)

    result = {"success": 1, 'n_obs': n_obs, 'n_genes': n_genes, 'top_genes': top_genes}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
