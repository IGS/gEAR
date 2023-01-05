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

import pandas as pd

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)
    n_genes = int(form.getvalue('n_genes'))
    compute_marker_genes = form.getvalue('compute_marker_genes')

    ana = geardb.Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id,
                          session_id=session_id, user_id=user.id)

    adata = ana.get_adata()

    dest_datafile_path = None
    cluster_method = 'louvain'

    if ana.type == 'primary':
        ana.type = 'user_unsaved'
        dest_datafile_path = ana.dataset_path()

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
    else:
        dest_datafile_path = ana.dataset_path()

    dest_datafile_dir = os.path.dirname(dest_datafile_path)

    if not os.path.exists(dest_datafile_dir):
        os.makedirs(dest_datafile_dir)

    # Previous steps have copied adata to raw, but if we got here via primary
    #  This probably hasn't happened.  Do it now.
    if not adata.raw:
        adata.raw = adata

    # Compute marker genes
    if compute_marker_genes == 'true':

        # TODO: add support for 'method' argument here, logreg, t-test, etc.'
        #print("DEBUG: For dataset {0} plotting with cluster_method: {1}".format(dataset_id, cluster_method), file=sys.stderr)
        sc.tl.rank_genes_groups(adata, cluster_method)
        adata.write(dest_datafile_path)

        ## I don't see how to get the save options to specify a directory
        # sc.settings.figdir = 'whateverpathyoulike' # scanpy issue #73
        os.chdir(dest_datafile_dir)

        #print("DEBUG: dataset_id:{0} wrote image file to: {1}".format(dataset_id, dest_datafile_dir), file=sys.stderr)

        # The sharey parameter here controls whether all axes have the same scale
        sc.pl.rank_genes_groups(adata, n_genes=n_genes, gene_symbols='gene_symbol', sharey=False, save='.png')

    df_json_str = pd.DataFrame(
        adata.uns['rank_genes_groups']['names']
    ).iloc[:n_genes].to_json(orient='split')

    """
    This section is to transform the JSON like this:

                {
                   columns: [
                      {label: 0},
                      {label: 1},
                      ...
                   ],
                   rows: [
                     {rowid: 0, columns: [ {label: 'Tagln'},
                                           {label: 'Arhgdib'}, ...
                                         ] },
                     {rowid: 1, columns: [ {label: 'Sparcl1'},
                                           {label: 'Hnrnph1'}, ...
                                         ] },
                   ]
                }
    """

    df_json = json.loads(df_json_str)
    tbl_json = {'columns': [], 'rows': []}

    for col in df_json['columns']:
        tbl_json['columns'].append({'label': col})

    row_idx = 0
    for row in df_json['data']:
        tbl_json['rows'].append({'rowid': df_json['index'][row_idx], 'columns':[]})

        for label in df_json['data'][row_idx]:
            gs = adata.raw.var.loc[label]
            if not isinstance(gs, str):
                gs = adata.raw.var.loc[label].get('gene_symbol')

            tbl_json['rows'][-1]['columns'].append({'label': gs})

        row_idx += 1

    # tbl_json_str = json.dumps(tbl_json)

    # this is for the gene group labels
    group_labels = []
    col_idx = 0
    for col_label in df_json['columns']:
        example_gene_str = df_json['data'][0][col_idx]

        gs = adata.raw.var.loc[example_gene_str]
        if not isinstance(gs, str):
            gs = adata.raw.var.loc[example_gene_str].get('gene_symbol')

        num_cells = adata.obs[adata.obs[cluster_method] == col_label][cluster_method].count()

        group_labels.append({'group_label':col_label, 'num_cells':num_cells, 'genes': gs})
        col_idx += 1

    result = {'success': 1, 'table': tbl_json, 'group_labels': group_labels, 'cluster_label': cluster_method}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
