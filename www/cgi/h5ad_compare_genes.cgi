#!/opt/bin/python3

"""

"""

import cgi
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.analysis import get_analysis, test_analysis_for_zarr

# this is needed so that we don't get TclError failures in the underlying modules
matplotlib.use('Agg')
sc.settings.verbosity = 0


def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    query_cluster = form.getvalue('query_cluster')
    reference_cluster = form.getvalue('reference_cluster')
    n_genes = int(form.getvalue('n_genes'))
    method = form.getvalue('method')
    corr_method = form.getvalue('corr_method')

    if form.getvalue('group_labels'):
        group_labels = json.loads(form.getvalue('group_labels'))
    else:
        group_labels = None

    ## correction method isn't valid for logistic regression
    if method == 'logreg':
        corr_method = None

    result = {"success": 0, "cluster_label": "", "table_json_f": None, "table_json_r": None }

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

    # Recent upgrade of scanpy/anndata/pandas modules have issues where the sc.pl.rank_genes_groups_violin function
    # fails as pandas throws an error saying the data is not 1-dimensional.  This only happens if the AnnData object is a dense matrix
    # My workaround is to force it to be sparse.
    kwargs = {}
    if not test_analysis_for_zarr(ana):
        kwargs["force_sparse"] = True

    # Try to get the adata in sparse format, if that fails, fall back to non-sparse
    try:
        adata = ana.get_adata(**kwargs)
    except Exception as e:
        # This is non-fatal
        print(f"INFO: {str(e)}. Switching to non-sparse representation.", file=sys.stderr)
        adata = ana.get_adata()

    cluster_method = 'louvain'

    if ana.type == 'primary':
        ## dirty hack for BICCN dataset customization
        # BICCN Mini-Atlas (Integrated)
        if dataset_id in ['4fbd43e2-f301-42d8-a61d-a5dd5bf720e7',
                          '7feba96b-816f-4154-8a3d-9f578acdf6c7',
                          'f70f6499-a319-4346-aca5-79d88e48744e',
                          '9682c015-9bb7-4307-bd7f-05abcd5d82b4']:
            cluster_method = 'joint_cluster_round4_annot'
                            # BICCN Mini-Atlas (Transcript)
        elif dataset_id in ['495d68ba-9c05-4f5a-8b39-fb62291ae205',
                            'bf902ed4-b163-4495-85fa-eea46e4d5670',
                            '0de50971-6518-4bee-a860-fc12cac8e623',
                            '72b8e735-6ac2-4e0d-919a-c989994f8be1',
                            '8eb9cd46-a5c8-40e6-983b-8809eb1201cc',
                            # Motor Cortex Cross Species
                            'aab7a4ac-a569-4d6e-9173-23a48950231f',
                            '4dffa5fd-e369-407c-abd1-d18b9f99a039',
                            'c4435959-a774-4689-865e-5ed0678aacdd',
                            '12545384-48b1-454b-a210-8781ca3e2230',
                            'dd7e9b4b-b127-4447-8d41-9fd16c954119',
                            'fcd8186e-f26f-493e-acc2-19e49edae191',
                            # Motor Cortex Merged Species
                            '071f042d-1b2c-4c2f-b330-248c9f5b58a1',
                            '2d7d38a6-7e17-40a7-bbc1-016afd2e5ec9',
                            'debbff92-dbe4-4b61-8cc8-b19d45ddf1d4'
        ]:
            cluster_method = 'subclass_label'

    # primary or public analyses won't be after this
    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path

    if group_labels:
        adata.rename_categories(cluster_method, group_labels)

    if reference_cluster == 'all-reference-clusters':
        if method == "logreg":
            error_msg = "Logistic regression method requires a specific reference cluster."
            print(error_msg, file=sys.stderr)
            result['success'] = 0
            result['error'] = error_msg
            sys.stdout = original_stdout
            print('Content-Type: application/json\n\n')
            print(json.dumps(result))
            return

        if corr_method is None:
            sc.tl.rank_genes_groups(adata, cluster_method, groups=[query_cluster], method=method)
        else:
            sc.tl.rank_genes_groups(adata, cluster_method, groups=[query_cluster], method=method, corr_method=corr_method)
    else:
        if corr_method is None:
            sc.tl.rank_genes_groups(adata, cluster_method, groups=[query_cluster], reference=reference_cluster, method=method)
        else:
            sc.tl.rank_genes_groups(adata, cluster_method, groups=[query_cluster], reference=reference_cluster, method=method, corr_method=corr_method)

    ## I don't see how to get the save options to specify a directory
    # sc.settings.figdir = 'whateverpathyoulike' # scanpy issue #73
    os.chdir(os.path.dirname(dest_datafile_path))

    adata.X = sparse.csc_matrix(adata.X)

    """
    One of the scanpy versions introduced a bug that was recently fixed, where pl.rank_genes_groups works
    but not pl.rank_genes_groups_<plot>.  I believe it is because Ensembl IDs are stored as gene_symbols for tl.rank_genes_groups
    and pl.rank_genes_groups looks for these ensembl_ids but pl.rank_genes_groups_<plot> is erroneously looking for the gene symbols instead

    Seems to work in scanpy 1.7.2 but is broke in 1.8.2 and was fixed in https://github.com/scverse/scanpy/pull/1529
    """

    ax = sc.pl.rank_genes_groups(adata, groups=[query_cluster],
                                 gene_symbols='gene_symbol', n_genes=n_genes, save="_comp_ranked.png")

    try:
        # Try 1.7.2 way first
        ax = sc.pl.rank_genes_groups_violin(adata, groups=query_cluster, use_raw = False,
                                            gene_symbols="gene_symbol", n_genes=n_genes, save="_comp_violin.png")
    except Exception:
        # Use gene names if that doesn't work
        gene_names = adata.var.loc[adata.uns["rank_genes_groups"]['names'][query_cluster]]["gene_symbol"][:n_genes].tolist()
        ax = sc.pl.rank_genes_groups_violin(adata, groups=query_cluster, use_raw = False,
                                            gene_symbols="gene_symbol", gene_names=gene_names, save="_comp_violin.png")

    result["success"] = 1
    result["cluster_label"] = cluster_method

    if not method == 'logreg':
        # Yuk.
        ensembl_id_list = np.concatenate(adata.uns['rank_genes_groups']['names'].tolist()).ravel().tolist()

        result['table_json_f'] = pd.DataFrame({
            'names': adata.var.loc[ensembl_id_list]['gene_symbol'].tolist(),
            'log2FC':adata.uns['rank_genes_groups']['logfoldchanges'],
            'P-value':adata.uns['rank_genes_groups']['pvals'],
            'FDR':adata.uns['rank_genes_groups']['pvals_adj']}).loc[:n_genes].to_json(orient='split')

    if reference_cluster != 'all-reference-clusters':
        if corr_method is None:
            sc.tl.rank_genes_groups(adata, cluster_method, groups=[reference_cluster], reference=query_cluster, method=method)
        else:
            sc.tl.rank_genes_groups(adata, cluster_method, groups=[reference_cluster], reference=query_cluster, method=method, corr_method=corr_method)

        ax = sc.pl.rank_genes_groups(adata, groups=[reference_cluster],
                                     gene_symbols='gene_symbol', n_genes=n_genes, save="_comp_ranked_rev.png")

        try:
            ax = sc.pl.rank_genes_groups_violin(adata, groups=reference_cluster, use_raw = False,
                                                gene_symbols="gene_symbol", n_genes=n_genes, save="_comp_violin_rev.png")
        except Exception:
            gene_names = adata.var.loc[adata.uns["rank_genes_groups"]['names'][reference_cluster]]["gene_symbol"][:n_genes].tolist()

            ax = sc.pl.rank_genes_groups_violin(adata, groups=reference_cluster, use_raw = False,
                                                gene_symbols="gene_symbol", gene_names=gene_names, save="_comp_violin_rev.png")

    if not method == 'logreg':
        # Yuk.
        ensembl_id_list = np.concatenate(adata.uns['rank_genes_groups']['names'].tolist()).ravel().tolist()

        # Copilot-proposed solution
        #import itertools
        #ensembl_id_list = list(itertools.chain.from_iterable(adata.uns['rank_genes_groups']['names'].tolist()))

        result['table_json_r'] = pd.DataFrame({
            'names': adata.var.loc[ensembl_id_list]['gene_symbol'].tolist(),
            'log2FC':adata.uns['rank_genes_groups']['logfoldchanges'],
            'P-value':adata.uns['rank_genes_groups']['pvals'],
            'FDR':adata.uns['rank_genes_groups']['pvals_adj']}).loc[:n_genes].to_json(orient='split')

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
