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

    plot_tsne = int(form.getvalue('plot_tsne'))
    plot_umap = int(form.getvalue('plot_umap'))

    ana = geardb.Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id,
                          session_id=session_id, user_id=user.id)

    resolution = float(form.getvalue('resolution'))
    compute_louvain = form.getvalue('compute_louvain')
    # Possible fix to tackle the error,
    # ValueError: Categorical categories must be unique
    # when renaming categories to these group labels
    # group_labels = list(set(json.loads(form.getvalue('group_labels'))))
    group_labels = json.loads(form.getvalue('group_labels'))

    adata = ana.get_adata()

    # primary or public analyses won't be after this
    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path()

    # Compute tSNE and plot
    if compute_louvain == 'true':
        # SAdkins note - Occasionally I run out of memory computing this step on Docker,
        # especially if I want to do downstream stuff.
        # If this happens, set 'flavor="igraph"' which uses a different package.
        sc.tl.louvain(adata, resolution=resolution)
        adata.write(dest_datafile_path)

    ## I don't see how to get the save options to specify a directory
    # sc.settings.figdir = 'whateverpathyoulike' # scanpy issue #73
    os.chdir(os.path.dirname(dest_datafile_path))

    # doing it like this puts the labels in the legend
    if len(group_labels) > 0:
        adata.rename_categories('louvain', group_labels)
        adata.obs['louvain'] = adata.obs['louvain'].astype(str)
        adata.write(dest_datafile_path)

        if plot_tsne == 1:
            ax = sc.pl.tsne(adata, color='louvain', legend_loc='on data', save="_louvain.png")
        if plot_umap == 1:
            ax = sc.pl.umap(adata, color='louvain', legend_loc='on data', save="_louvain.png")
    else:
        if plot_tsne == 1:
            ax = sc.pl.tsne(adata, color='louvain', save="_louvain.png")
        if plot_umap == 1:
            ax = sc.pl.umap(adata, color='louvain', save="_louvain.png")

    result = {'success': 1}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

