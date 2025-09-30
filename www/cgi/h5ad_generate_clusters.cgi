#!/opt/bin/python3

"""

"""

import cgi
import json
import os
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

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')

    plot_tsne = int(form.getvalue('plot_tsne'))
    plot_umap = int(form.getvalue('plot_umap'))

    result = {"success": 0, "group_labels":""}

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
    resolution = float(form.getvalue('resolution'))
    compute_clusters = form.getvalue('compute_clusters')
    cluster_info = json.loads(form.getvalue("cluster_info"))    # "old_label", "new_label", "keep"

    adata = ana.get_adata()

    # primary or public analysis should not be overwritten
    # this will alter the analysis object save destination
    if ana.type == 'primary' or ana.type == 'public':
        ana.type = 'user_unsaved'

    dest_datafile_path = ana.dataset_path

    # Compute tSNE and plot
    if compute_clusters == 'true':

        adata.obs.drop(columns=["louvain", "orig_louvain"], errors="ignore", inplace=True)

        try:
            # Added flavor and n_iterations to address warnings about future defaults
            sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2)

            # rename the leiden column to louvain to not break things elsewhere
            # ? perhaps we should rename as "clustering" or something more generic
            adata.obs.rename(columns={"leiden":"louvain"}, inplace=True)
        except Exception as e:
            print("Switching to louvain algorithm, leiden failed", file=sys.stderr)
            print("Error: ", e, file=sys.stderr)
            sc.tl.louvain(adata, resolution=resolution, flavor="igraph")

        adata.obs["orig_louvain"] = adata.obs["louvain"].astype(int)   # Copy cluster ID so it's easier to rename categories
        adata.write(dest_datafile_path)

    ## I don't see how to get the save options to specify a directory
    # sc.settings.figdir = 'whateverpathyoulike' # scanpy issue #73
    os.chdir(os.path.dirname(dest_datafile_path))

    group_labels = []

    # doing it like this puts the labels in the legend
    if len(cluster_info) > 0:

        # If this is an older louvain analysis, make this mapping column if it does not exist
        if "orig_louvain" not in adata.obs:
            old_label2index = dict()
            for idx, cluster in enumerate(cluster_info):
                old_label2index[cluster["old_label"]] = idx
            adata.obs["orig_louvain"] = adata.obs["louvain"].map(old_label2index)

        # ? I think mapping old label to new label w/o using index would have worked fine, but this cleans up the logic for me
        # Temporarily make the louvain IDs the index numbers
        adata.obs["louvain"] = adata.obs["orig_louvain"]

        # Filter only the clusters we want to use
        kept_indexes = list(filter(lambda i: cluster_info[i]["keep"], range(len(cluster_info))))
        adata = adata[adata.obs["louvain"].isin(kept_indexes), :]

        # Create mapping of original cluster IDs and new labels. Clusters will merge on duplicated labels
        idx2new_label = dict()
        for idx, cluster in enumerate(cluster_info):
            idx2new_label[idx] = cluster["new_label"]
        adata.obs["louvain"] = adata.obs["louvain"].map(idx2new_label)

        # Create new cluster IDs and labels. Assumes that running this script again will preserve the order
        # Duplicating "h5ad_find_marker_genes" group labels structure here, so we can re-render the html table
        label2idx = dict()
        deduped_group_labels = adata.obs["louvain"].unique().tolist()
        for idx, label in enumerate(deduped_group_labels):
            num_cells = adata.obs[adata.obs["louvain"] == label]["louvain"].count()
            if not num_cells:
                continue
            group_labels.append({'group_label':idx, 'num_cells':num_cells, 'genes': label})
            label2idx[label] = idx

        # Ensure orig_louvain is parallel to the group_labels, so relabeling uses the correct cluster numbers
        if not (len(cluster_info) == len(deduped_group_labels)):
            adata.obs["orig_louvain"] = adata.obs["louvain"].map(label2idx)

        # Ensure louvain is string for downstream uses
        adata.obs['louvain'] = adata.obs['louvain'].astype(str)

        adata.write(dest_datafile_path)


        # Rename "louvain" to a generic "clustering" for consistency
        adata.obs.rename(columns={"louvain":"clustering"}, inplace=True)

        if plot_tsne == 1:
            ax = sc.pl.tsne(adata, color='clustering', legend_loc='on data', save="_clustering.png")
        if plot_umap == 1:
            ax = sc.pl.umap(adata, color='clustering', legend_loc='on data', save="_clustering.png")
    else:
        adata.obs.rename(columns={"louvain":"clustering"}, inplace=True)
        if plot_tsne == 1:
            ax = sc.pl.tsne(adata, color='clustering', save="_clustering.png")
        if plot_umap == 1:
            ax = sc.pl.umap(adata, color='clustering', save="_clustering.png")

    result["success"] = 1
    result["group_labels"] = group_labels

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
