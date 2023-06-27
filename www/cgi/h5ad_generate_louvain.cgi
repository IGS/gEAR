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

import pandas as pd

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
    group_labels = json.loads(form.getvalue('group_labels'))

    adata = ana.get_adata()

    # primary or public analyses won't be after this
    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path()

    # Compute tSNE and plot
    if compute_louvain == 'true':
        # NOTE - Occasionally I run out of memory computing this step on Docker,
        # especially if I want to do downstream stuff.
        # If this happens, set 'flavor="igraph"' which uses a different package.
        sc.tl.louvain(adata, resolution=resolution, flavor="igraph")
        adata.obs["orig_louvain"] = adata.obs["louvain"]   # Copy order so it's easier to rename categories
        adata.write(dest_datafile_path)

    ## I don't see how to get the save options to specify a directory
    # sc.settings.figdir = 'whateverpathyoulike' # scanpy issue #73
    os.chdir(os.path.dirname(dest_datafile_path))

    # doing it like this puts the labels in the legend
    if len(group_labels) > 0:

        # NOTE: This is not backwards compatible with louvain computations before this was added.
        # For those, renaming labels 2+ times requires a full analyses rerun (to reset louvain)
        if "orig_louvain" in adata.obs:
            adata.obs["louvain"] = adata.obs["orig_louvain"]

        # Create mapping of original cluster IDs and new labels. Clusters will merge on duplicated labels
        idx_label_map = dict()
        for idx, label in enumerate(group_labels):
            str_idx = str(idx)  # sc.tl.louvain always saves clusters as strings
            idx_label_map[str_idx] = label
        adata.obs["louvain"] = adata.obs["louvain"].map(idx_label_map)

        # Create new cluster IDs and labels. Assumes that running this script again will preserve the order
        # Duplicating "h5ad_find_marker_genes" group labels structure here, so we can re-render the html table
        new_group_labels = []
        label2idx = dict()
        deduped_group_labels = list(set(group_labels))
        for idx, label in enumerate(deduped_group_labels):
            num_cells = adata.obs[adata.obs["louvain"] == label]["louvain"].count()
            new_group_labels.append({'group_label':idx, 'num_cells':num_cells, 'genes': label})
            label2idx[label] = str(idx)
            # ? Can we eliminate making as string since and use ints in "louvain" and "orig_louvain"

        # Ensure orig_louvain is parallel to the group_labels, so relabeling uses the correct cluster numbers
        if not len(group_labels) == len(deduped_group_labels):
            adata.obs["orig_louvain"] = adata.obs["louvain"].map(label2idx)
        group_labels = new_group_labels

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

    result = {'success': 1, "group_labels":group_labels}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()

