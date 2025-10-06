#!/opt/bin/python3

"""
This script is a CGI handler for adding or updating primary analysis metadata to a dataset in the gEAR platform.
It performs the following main tasks:
- Authenticates the user via session ID.
- Locates the dataset using a share UID.
- Loads the primary analysis object and its associated settings JSON.
- Detects the presence of tSNE, UMAP, and clustering results in the dataset's AnnData object.
- Updates the AnnData object and analysis JSON to ensure tSNE, UMAP, and clustering results are present and properly recorded.
- Writes changes to the AnnData file and analysis JSON if modifications are made.
- Returns a JSON response indicating success or failure.

Key Functions:
- main(): Orchestrates the process of updating the primary analysis for a dataset.
- add_clustering_analysis(adata): Adds clustering information to the AnnData object if detected.
- add_tsne_analysis(adata): Adds tSNE coordinates to the AnnData object if detected.
- add_umap_analysis(adata): Adds UMAP coordinates to the AnnData object if detected.
- detect_clustering(adata): Checks if clustering columns are present in AnnData.obs.
- detect_tsne(adata): Checks if tSNE coordinate columns are present in AnnData.obs.
- detect_umap(adata): Checks if UMAP coordinate columns are present in AnnData.obs.
- has_clustering(adata): Checks if clustering results are already present in AnnData.obs.
- has_tsne(adata): Checks if tSNE results are already present in AnnData.obsm.
- has_umap(adata): Checks if UMAP results are already present in AnnData.obsm.

Constants:
- DATASET_BASE_DIR: Base directory for datasets.
- VALID_CLUSTER_COLUMN_NAMES: List of valid column names for clustering.
- VALID_TSNE_PAIRS: List of valid tSNE coordinate column pairs.
- VALID_UMAP_PAIRS: List of valid UMAP coordinate column pairs.

The script is intended to be run as a CGI script during the upload process and outputs a JSON response.
"""


import cgi
import json
import os
import sys
import typing
from pathlib import Path

import scanpy as sc

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

gear_root_path = Path(__file__).resolve().parents[2]

lib_path = gear_root_path.joinpath('lib')
sys.path.insert(0, str(lib_path))

import geardb
from gear.analysis import H5adAdapter, ZarrAdapter

sc.settings.verbosity = 0

if typing.TYPE_CHECKING:
    # This allows type-checkers to resolve types without importing the actual modules at runtime.
    # To avoid having runtime errors, enclose the typing in quotes (AKA forward-reference)
    from anndata import AnnData

DATASET_BASE_DIR = '{}/www/datasets'.format(gear_root_path)
VALID_CLUSTER_COLUMN_NAMES = ['cluster', 'cell_type', 'cluster_label', 'subclass_label', 'joint_cluster_round4_annot']
VALID_TSNE_PAIRS = [['tSNE_1', 'tSNE_2'], ['tSNE1', 'tSNE2'], ['tsne1_combined', 'tsne2_combined'],
                    # These are all carlo's custom ones.  Need to resolve this a different way later
                    ['PC1%var6.14', 'PC2%var1.79']]
VALID_UMAP_PAIRS = [['uMAP_1', 'uMAP_2'], ['uMAP1', 'uMAP2'],
                    ['UMAP_1', 'UMAP_2'], ['UMAP1', 'UMAP2']]


user_upload_file_base = gear_root_path / 'www' / 'uploads' / 'files'

def main() -> dict:
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    share_uid = form.getvalue('share_uid')
    dataset_format = form.getvalue('dataset_format', 'single-cell-rnaseq')

    result = {"success": 0, "message": "", 'perform_primary_analysis': True}

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        print(result['message'], file=sys.stderr)
        return result

    if not share_uid:
        result['message'] = 'No share ID found for dataset. Please try upload again.'
        print(result['message'], file=sys.stderr)
        return result

    print("Processing share ID: {0}".format(share_uid))

    # Load the metadata, to get spatial information
    dataset_upload_dir = user_upload_file_base / session_id / share_uid
    metadata_file = dataset_upload_dir / "metadata.json"
    if not metadata_file.is_file():
        result['message'] = 'Metadata file not found.'
        print(result['message'], file=sys.stderr)
        return result

    with open(metadata_file, 'r') as f:
        metadata = json.load(f)
        dataset_id = metadata.get('dataset_uid', '')

    # If the dataset type is not single-cell-rnaseq or spatial, exit gracefully
    if metadata.get("dataset_type") not in ['single-cell-rnaseq', 'spatial']:
        result['success'] = 1
        result['perform_primary_analysis'] = False
        result['message'] = 'This dataset format cannot be run in the single-cell workbench. Exiting gracefully.'
        print(result['message'], file=sys.stderr)
        return result

    # Update that primary analysis can be performed on this dataset
    metadata["perform_primary_analysis"] = True
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=4)

    # Load the analysis JSON or create from template
    analysis_json_path = dataset_upload_dir / "analysis_pipeline.json"
    if analysis_json_path.is_file():
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
    analysis_json["dataset"]["id"] = dataset_id

    # Figure out how to retrieve the AnnData object based on the file type
    kwargs = {}
    if dataset_format == "spatial":
        upload_file = dataset_upload_dir / f"{share_uid}.zarr"
        if not upload_file.is_dir():
            result['message'] = f"Spatial dataset file not found for share ID: {share_uid}"
            print(result['message'], file=sys.stderr)
            return result
        adapter = ZarrAdapter(upload_file)
    else:
        upload_file = dataset_upload_dir / f"{share_uid}.h5ad"
        if not upload_file.is_file():
            result['message'] = f"Dataset file not found for share ID: {share_uid}"
            print(result['message'], file=sys.stderr)
            return result
        adapter = H5adAdapter(upload_file)
        kwargs["backed"] = True

    adata = adapter.get_adata(**kwargs)

    h5ad_changes_made = False
    json_changes_made = False
    # Does the metadata file exist?  If not, we need to create it
    if not metadata_file.is_file():
        json_changes_made = True

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
