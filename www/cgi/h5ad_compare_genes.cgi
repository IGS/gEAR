#!/opt/bin/python3

"""

"""

import cgi, json
import os, sys, re
import numpy as np

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

from scipy import sparse

import pandas as pd

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)
    query_cluster = form.getvalue('query_cluster')
    reference_cluster = form.getvalue('reference_cluster')
    n_genes = int(form.getvalue('n_genes'))
    method = form.getvalue('method')
    corr_method = form.getvalue('corr_method')
    compute_gene_comparison = int(form.getvalue('compute_gene_comparison'))

    if form.getvalue('group_labels'):
        group_labels = json.loads(form.getvalue('group_labels'))
    else:
        group_labels = None

    ## correction method isn't valid for logistic regression
    if method == 'logreg':
        corr_method = None

    ana = geardb.Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id,
                          session_id=session_id, user_id=user.id)

    # Recent upgrade of scanpy/anndata/pandas modules have issues where the sc.pl.rank_genes_groups_violin function
    # fails as pandas throws an error saying the data is not 1-dimensional.  This only happens if the AnnData object is a dense matrix
    # My workaround is to force it to be sparse.
    adata = ana.get_adata(force_sparse=True)
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

    dest_datafile_path = ana.dataset_path()

    if group_labels:
        adata.rename_categories(cluster_method, group_labels)

    if reference_cluster == 'all-reference-clusters':
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
    
    ax = sc.pl.rank_genes_groups(adata, groups=[query_cluster],
                                 gene_symbols='gene_symbol', n_genes=n_genes, save="_comp_ranked.png")
    ax = sc.pl.rank_genes_groups_violin(adata, groups=query_cluster, use_raw=False,
                                        gene_symbols='gene_symbol', n_genes=n_genes, save="_comp_violin.png")

    result = {'success': 1, 'cluster_label': cluster_method}

    if method == 'logreg':
        result['table_json_f'] = ''
    else:
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
        ax = sc.pl.rank_genes_groups_violin(adata, groups=reference_cluster, use_raw = False,
                                            gene_symbols='gene_symbol', n_genes=n_genes, save="_comp_violin_rev.png")

    if method == 'logreg':
        result['table_json_r'] = ''
    else:
        # Yuk.
        ensembl_id_list = np.concatenate(adata.uns['rank_genes_groups']['names'].tolist()).ravel().tolist()

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

