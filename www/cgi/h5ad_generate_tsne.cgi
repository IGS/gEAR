#!/opt/bin/python3

"""

"""

import cgi
import json
import os
import re
import sys

import matplotlib
import scanpy as sc

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.analysis import get_analysis

# this is needed so that we don't get TclError failures in the underlying modules

matplotlib.use('Agg')

sc.settings.verbosity = 0

def normalize_genes_to_color(gene_list, chosen_genes):
    """Convert to case-insensitive.  Also will not add chosen gene if not in gene list."""
    case_insensitive_genes = [g for cg in chosen_genes for g in gene_list if cg.lower() == g.lower()]
    return case_insensitive_genes


def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')

    n_pcs = int(form.getvalue('n_pcs'))
    n_neighbors = int(form.getvalue('n_neighbors'))
    random_state = int(form.getvalue('random_state'))
    genes_to_color = form.getvalue('genes_to_color')
    use_scaled = form.getvalue('use_scaled')

    compute_neighbors = int(form.getvalue('compute_neighbors'))
    compute_tsne = int(form.getvalue('compute_tsne'))
    compute_umap = int(form.getvalue('compute_umap'))

    plot_tsne = int(form.getvalue('plot_tsne'))
    plot_umap = int(form.getvalue('plot_umap'))

    result = {"success": 0}

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        print("No dataset found with that ID.", file=sys.stderr)
        result['success'] = 0
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return
    is_spatial = ds.dtype == "spatial"

    analysis_obj = None
    if analysis_id or analysis_type:
        analysis_obj = {
            'id': analysis_id if analysis_id else None,
            'type': analysis_type if analysis_type else None,
        }

    try:
        ana = get_analysis(analysis_obj, dataset_id, session_id, is_spatial=is_spatial)
    except Exception:
        print("Analysis for this dataset is unavailable.", file=sys.stderr)
        result['success'] = 0
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        return

    if genes_to_color:
        genes_to_color = genes_to_color.replace(' ', '')

        if ',' in genes_to_color:
            genes_to_color = genes_to_color.split(',')
        else:
            genes_to_color = [genes_to_color]

    adata = ana.get_adata()

    # primary or public analysis should not be overwritten
    # this will alter the analysis object save destination
    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path

    if compute_neighbors == 1:
        sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    if compute_tsne == 1:
        sc.tl.tsne(adata, n_pcs=n_pcs, random_state=random_state)

    if compute_umap == 1:
        sc.tl.umap(adata, maxiter=500)

    # If any of the above steps were done, save the adata object
    if compute_neighbors == 1 or compute_tsne == 1 or compute_umap == 1:
        adata.write(dest_datafile_path)
    else:
        # Get from the dest_datafile_path
        adata = ana.get_adata()

    ## I don't see how to get the save options to specify a directory
    # sc.settings.figdir = 'whateverpathyoulike' # scanpy issue #73
    os.chdir(os.path.dirname(dest_datafile_path))

    missing_gene = None

    if genes_to_color:
        # Catch the error if any gene names are passed which aren't in the dataset
        try:
            adata.var = adata.var.reset_index().set_index('gene_symbol')
            # We also need to change the adata's Raw var dataframe
            # We can't explicitly reset its index so we reinitialize it with
            # the newer adata object.
            # https://github.com/theislab/anndata/blob/master/anndata/base.py#L1020-L1022
            if adata.raw is not None:
                adata.raw = adata

            gene_symbols = adata.var.index.tolist()
            genes_to_color = normalize_genes_to_color(gene_symbols, genes_to_color)

            # This can error like: ValueError: key "RFX7" is invalid! specify valid sample annotation
            # original color map: RdBu_r
            if use_scaled == 'true':
                if plot_tsne == 1:
                    sc.pl.tsne(adata, color=genes_to_color, color_map='YlOrRd', use_raw=False, save=".png")

                if plot_umap == 1:
                    sc.pl.umap(adata, color=genes_to_color, color_map='YlOrRd', use_raw=False, save=".png")
            else:
                if plot_tsne == 1:
                    sc.pl.tsne(adata, color=genes_to_color, color_map='YlOrRd', save=".png")

                if plot_umap == 1:
                    sc.pl.umap(adata, color=genes_to_color, color_map='YlOrRd', save=".png")
        except ValueError as err:
            # scanpy seems to change this error string every release
            #print("DEBUG: error string:{0}".format(str(err)), file=sys.stderr)
            # DEBUG: error string:Given 'color': foobar is not a valid observation or var. Valid observations are: Index(['n_genes', 'n_counts'], dtype='object')
            m = re.search("\: (.+?) is not a valid", str(err))
            if m:
                missing_gene = m.groups(1)
            else:
                missing_gene = 'Unknown'
    else:
        if plot_tsne == 1:
            sc.pl.tsne(adata, color_map='YlOrRd', save=".png")

        if plot_umap == 1:
            sc.pl.umap(adata, color_map='YlOrRd', save=".png")

    if missing_gene is None:
        result = {'success': 1, 'missing_gene': ''}
    else:
        result = {'success': 0, 'missing_gene': missing_gene}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

