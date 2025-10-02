#!/opt/bin/python3

"""

Many datasets exist in the system already with columnar elements which
make up a user-generated analysis, such as dimensionality reduction (tSNE)
and clustering.

This script reads these and transforms them into a structure within HDF5
which mimics what scanpy produces.  This allows them to be used throughout
the rest of the interface.

Autodetection isn't perfect.  Currently it checks for these columns for
clustering identification:

VALID_CLUSTER_COLUMN_NAMES = ['cluster', 'cell_type', 'cluster_label', 'subclass_label']

And these pairs for tSNE:

VALID_TSNE_PAIRS = [['tSNE_1', 'tSNE_2'], ['tSNE1', 'tSNE2'], ['tsne1_combined', 'tsne2_combined']]
VALID_UMAP_PAIRS = [['uMAP_1', 'uMAP_2'], ['uMAP1', 'uMAP2'],
                    ['UMAP_1', 'UMAP_2'], ['UMAP1', 'UMAP2']]


"""

import cgi
import json
import os
import re
import sys
import typing
from pathlib import Path

import scanpy as sc

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = Path(__file__).resolve().parents[2].joinpath('lib')
sys.path.insert(0, str(lib_path))

import geardb
from gear.analysis import get_primary_analysis, test_analysis_for_zarr

sc.settings.verbosity = 0

if typing.TYPE_CHECKING:
    # This allows type-checkers to resolve types without importing the actual modules at runtime.
    # To avoid having runtime errors, enclose the typing in quotes (AKA forward-reference)
    from anndata import AnnData
    from gear.analysis import Analysis, SpatialAnalysis

# These should match in cgi/get_embedded_tsne_display.cgi
gear_root_path = lib_path.parents[1]

DATASET_BASE_DIR = '{}/www/datasets'.format(gear_root_path)
VALID_CLUSTER_COLUMN_NAMES = ['cluster', 'cell_type', 'cluster_label', 'subclass_label', 'joint_cluster_round4_annot']
VALID_TSNE_PAIRS = [['tSNE_1', 'tSNE_2'], ['tSNE1', 'tSNE2'], ['tsne1_combined', 'tsne2_combined'],
                    # These are all carlo's custom ones.  Need to resolve this a different way later
                    ['PC1%var6.14', 'PC2%var1.79']]
VALID_UMAP_PAIRS = [['uMAP_1', 'uMAP_2'], ['uMAP1', 'uMAP2'],
                    ['UMAP_1', 'UMAP_2'], ['UMAP1', 'UMAP2']]


def main() -> dict:
    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')

    result = {"success": 0}

    if not dataset_id:
        print("No dataset ID provided. Processing all datasets.", file=sys.stderr)
        return result


    print("Processing dataset ID: {0}".format(dataset_id))

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        print("No dataset found with that ID.", file=sys.stderr)
        result['success'] = 0
        return result

    is_spatial = ds.dtype == "spatial"

    try:
        ana = get_primary_analysis(dataset_id, is_spatial=is_spatial)
    except Exception:
        print("Analysis for this dataset is unavailable.", file=sys.stderr)
        result['success'] = 0
        return result

    ana.vetting = 'owner'

    analysis_json_path = ana.settings_path

    if os.path.exists(analysis_json_path):
        with open(analysis_json_path) as json_in:
            analysis_json = json.load(json_in)
    else:
        JSON_PATH = gear_root_path / "data" / "analysis_pipeline_template.json"
        with open(JSON_PATH) as json_in:
            analysis_json = json.load(json_in)

    # Analysis ID and dataset ID are the same for this type
    analysis_json["id"] = dataset_id
    analysis_json["type"] = "primary"
    analysis_json["label"] = "Primary analysis"
    analysis_json["dataset_id"] = dataset_id
    analysis_json["dataset"]["id"] = dataset_id

    # Load the AnnData object
    kwargs = {}
    if not test_analysis_for_zarr(ana):
        kwargs["backed"] = True

    adata = ana.get_adata(**kwargs)

    h5ad_changes_made = False
    json_changes_made = False

    tsne_detected = detect_tsne(adata)
    if tsne_detected:
        if 'tsne' not in analysis_json or not analysis_json['tsne']['tsne_calculated']:
            json_changes_made = True

        analysis_json['tsne']['tsne_calculated'] = True
        analysis_json['tsne']['plot_tsne'] = 1
        if not has_tsne(adata):
            print("\tAdding tSNE analysis")
            add_tsne_analysis(adata)
            h5ad_changes_made = True

    umap_detected = detect_umap(adata)
    if umap_detected:
        if 'tsne' not in analysis_json or 'tsne_calculated' not in analysis_json['tsne'] or not analysis_json['tsne']['tsne_calculated']:
            json_changes_made = True

        # the parent key 'tsne' here should really be renamed to 'dimred' or something like it
        analysis_json['tsne']['umap_calculated'] = True
        analysis_json['tsne']['plot_umap'] = 1
        if not has_umap(adata):
            print("\tAdding UMAP analysis")
            add_umap_analysis(adata)
            h5ad_changes_made = True

    clustering_detected = detect_clustering(adata)
    if clustering_detected:
        if 'clustering' not in analysis_json or not analysis_json['clustering']['calculated']:
            json_changes_made = True

        analysis_json['clustering']['calculated'] = True
        if not has_clustering(adata):
            print("\tAdding clustering analysis")
            add_clustering_analysis(adata)
            h5ad_changes_made = True

    if h5ad_changes_made:
        # 482cb81f-4816-e4e3-a3fe-514707b847d8.h5ad
        print("\tWriting new H5AD")
        adata.write()

    # Does the JSON need to be updated?
    if json_changes_made:
        print("\tWriting new JSON")
        with open(analysis_json_path, 'w') as outfile:
            json.dump(analysis_json, outfile, indent=3)

    result['success'] = 1
    return result

def add_clustering_analysis(adata: "AnnData") -> None:
    cols = adata.obs.columns.tolist()

    for vname in VALID_CLUSTER_COLUMN_NAMES:
        if vname in cols:
            user_defined_cluster_names = adata.obs[vname].astype('category')
            adata.obs['louvain'] = user_defined_cluster_names
            return

def add_tsne_analysis(adata: "AnnData") -> None:
    cols = adata.obs.columns.tolist()

    for pair in VALID_TSNE_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            adata.obsm['X_tsne'] = adata.obs[[pair[0], pair[1]]].values
            return

def add_umap_analysis(adata: "AnnData") -> None:
    cols = adata.obs.columns.tolist()

    for pair in VALID_UMAP_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            adata.obsm['X_umap'] = adata.obs[[pair[0], pair[1]]].values
            return

def detect_clustering(adata: "AnnData") -> bool:
    """
    Looks first for 'cluster', then 'cell_type'
    """
    cols = adata.obs.columns.tolist()

    for vname in VALID_CLUSTER_COLUMN_NAMES:
        if vname in cols:
            return True

    return False

def detect_tsne(adata: "AnnData") -> bool:
    """
    Looks for the combination of pairs in VALID_TSNE_PAIRS
    """
    cols = adata.obs.columns.tolist()

    for pair in VALID_TSNE_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            return True

    return False

def detect_umap(adata: "AnnData") -> bool:
    """
    Looks for the combination of pairs in VALID_UMAP_PAIRS
    """
    cols = adata.obs.columns.tolist()

    for pair in VALID_UMAP_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            return True

    return False

def has_clustering(adata: "AnnData") -> bool:
    cols = adata.obs.columns.tolist()
    # "louvain" is a legacy name when we used the scanpy
    # louvain analysis instead of the modern leiden one.
    if 'louvain' in cols:
        return True
    else:
        return False

def has_tsne(adata: "AnnData") -> bool:
    try:
        if type(adata.obsm['X_tsne']):
            return True
    except Exception:
        pass

    return False

def has_umap(adata: "AnnData") -> bool:
    try:
        if type(adata.obsm['X_umap']):
            return True
    except Exception:
        pass

    return False

if __name__ == '__main__':
    result = main()
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))
