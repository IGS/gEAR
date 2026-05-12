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

def create_broken_stick_model(num_pcs, sum_pcs):
    """Create standard broken stick model and return ratios."""
    bsm = [ 0 for i in range(0, num_pcs) ]
    bsm[0] = 1 / num_pcs
    for i in range(1,num_pcs):
        bsm[i] = bsm[i-1] + (1 / (num_pcs - i))
    bsm_ratio = [elem / num_pcs for elem in bsm]
    norm_bsm_ratio = [ sum_pcs * elem for elem in bsm_ratio]
    return norm_bsm_ratio[::-1]  # Reverse list to have largest ratio first

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

    genes_to_color = form.getvalue('genes_to_color')
    compute_pca = form.getvalue('compute_pca')

    if genes_to_color:
        genes_to_color = genes_to_color.rstrip().replace(' ', '')

    if genes_to_color and ',' in genes_to_color:
        genes_to_color = genes_to_color.split(',')
    elif genes_to_color:
        genes_to_color = [genes_to_color]
    else:
        genes_to_color = []

    adata = ana.get_adata()

    # primary or public analysis should not be overwritten
    # this will alter the analysis object save destination
    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path
    dest_directory = os.path.dirname(dest_datafile_path)

    if not os.path.exists(dest_directory):
        os.makedirs(dest_directory)

    # Compute PCA and make a scatter plot.
    if compute_pca == 'true':
        sc.tl.pca(adata, svd_solver='arpack')
        adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat R
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

            gene_symbols = adata.var.index.tolist()
            genes_to_color = normalize_genes_to_color(gene_symbols, genes_to_color)

            ax = sc.pl.pca(adata, color=genes_to_color, save=".png")
        except ValueError as err:
            m = re.search("\: (.+?) is not a valid", str(err))
            if m:
                missing_gene = m.groups(1)
            else:
                missing_gene = 'Unknown'
    else:
        ax = sc.pl.pca(adata, save=".png")

    if missing_gene is None:
        sc.pl.pca_variance_ratio(adata, log=True, save=".png")

        """ Broken stick model code.
        # Add variance ratio (and other info) to adata object
        sc.tl.pca(adata, n_comps=20)
        variance_ratios = adata.uns['pca']['variance_ratio'].tolist()
        num_pcs = len(variance_ratios)

        # The sum of the ratios will most likely not equal 1.
        # scanpy.tl.pca computes as many PCs as possible but only stores data for the first n_comps
        sum_pcs = sum(variance_ratios)

        bcm_ratios = create_broken_stick_model(num_pcs, sum_pcs)

        x = [i+1 for i in range(0, num_pcs)]
        pc_labels = [str(i) for i in x]

        fig, ax = plt.subplots()
        ax.bar(x, variance_ratios, label='PC Variance')
        ax.plot(x, bcm_ratios, 'r-o', label='Broken Stick Model')  # Red line plot w/ markers
        ax.set_xticks(x)
        ax.set_xticklabels(pc_labels)
        ax.set_xlabel("Principal Component")
        ax.set_ylabel("Variance Ratio")
        ax.legend(loc='upper right')
        fig.savefig('./figures/pca_variance_ratio.png', transparent=True, bbox_inches="tight")
        """

        result = {'success': 1, 'missing_gene': ''}
    else:
        result = {'success': 0, 'missing_gene': missing_gene}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

