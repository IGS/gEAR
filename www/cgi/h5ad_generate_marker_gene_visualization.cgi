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

def normalize_marker_genes(gene_list, chosen_genes):
    """Convert to case-insensitive.  Also will not add chosen gene if not in gene list."""
    case_insensitive_genes = [g for cg in chosen_genes for g in gene_list if cg.lower() == g.lower()]
    return case_insensitive_genes

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)

    ana = geardb.Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id,
                          session_id=session_id, user_id=user.id)

    marker_genes = json.loads(form.getvalue('marker_genes'))

    # client may send empty string when generating marker
    # gene visualizations more than once
    marker_genes = list(filter(lambda x: x != '', marker_genes))

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

    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path()
    dest_datafile_dir = os.path.dirname(dest_datafile_path)

    if not os.path.exists(dest_datafile_dir):
        os.makedirs(dest_datafile_dir)

    ## I don't see how to get the save options to specify a directory
    # sc.settings.figdir = 'whateverpathyoulike' # scanpy issue #73
    os.chdir(dest_datafile_dir)

    # Reset index to gene symbol to ensure gene names are the plot labels
    adata.var.reset_index(inplace=True)
    adata.var.set_index('gene_symbol', inplace=True)

    # Deduplicate gene_symbols
    adata = adata[:, adata.var.index.duplicated() == False]

    adata.var_names_make_unique()
    gene_symbols = adata.var.index.tolist()
    marker_genes = normalize_marker_genes(gene_symbols, marker_genes)

    sc.pl.dotplot(adata, marker_genes, groupby=cluster_method, use_raw=False, save='goi.png')
    sc.pl.stacked_violin(adata, marker_genes, groupby=cluster_method, use_raw=False, save='goi.png')

    result = {'success': 1}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

