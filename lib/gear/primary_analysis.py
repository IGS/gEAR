"""
Module centered around adding a primary analysis to a dataset.  This includes:
- Detecting whether tSNE, UMAP, or clustering analyses have already been performed based
    on the presence of certain columns in the AnnData object
- Applies to single-cell-rnaseq datasets and spatial datasets.
"""

import json
import shutil
import typing
from pathlib import Path

import scanpy as sc

from .analysis import H5adAdapter, ZarrAdapter

if typing.TYPE_CHECKING:
    # This allows type-checkers to resolve types without importing the actual modules at runtime.
    # To avoid having runtime errors, enclose the typing in quotes (AKA forward-reference)
    from anndata import AnnData

gear_root_path = Path(__file__).resolve().parents[2]

VALID_CLUSTER_COLUMN_NAMES = ['cluster', 'cell_type', 'cluster_label', 'subclass_label', 'joint_cluster_round4_annot']
VALID_TSNE_PAIRS = [['tSNE_1', 'tSNE_2'], ['tSNE1', 'tSNE2'], ['tsne1_combined', 'tsne2_combined'],
                    # These are all carlo's custom ones.  Need to resolve this a different way later
                    ['PC1%var6.14', 'PC2%var1.79']]
VALID_UMAP_PAIRS = [['uMAP_1', 'uMAP_2'], ['uMAP1', 'uMAP2'],
                    ['UMAP_1', 'UMAP_2'], ['UMAP1', 'UMAP2']]

sc.settings.verbosity = 0
sc.settings.autosave = True
sc.settings.figdir = '/tmp' # for the composition plots

class PrimaryAnalysisProcessingError(Exception):
    """Custom exception for errors during primary analysis processing."""
    pass

def add_primary_analysis_to_dataset(dataset_id, share_id, staging_dir, dataset_format):

    # Load the analysis JSON or create from template
    analysis_json_path = staging_dir / "analysis_pipeline.json"
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
        upload_file = staging_dir / f"{share_id}.zarr"
        if not upload_file.is_dir():
            raise PrimaryAnalysisProcessingError(f"Spatial dataset file not found for share ID: {share_id}")
        adapter = ZarrAdapter(upload_file)
    else:
        upload_file = staging_dir / f"{share_id}.h5ad"
        if not upload_file.is_file():
            raise PrimaryAnalysisProcessingError(f"Dataset file not found for share ID: {share_id}")
        adapter = H5adAdapter(upload_file)
        kwargs["backed"] = True

    adata = adapter.get_adata(**kwargs)

    # Create some initial composition plots
    create_composition_plots(adata, staging_dir, dataset_format == "spatial")

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
    return True

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

def create_composition_plots(adata: "AnnData", dataset_path: str, is_spatial: bool) -> None:
    # Create the pathnames for the images
    extension = ".h5ad"
    if is_spatial:
        extension = ".zarr"

    violin_image_path = str(dataset_path).replace(extension, '.prelim_violin.png')
    scatter_image_path = str(dataset_path).replace(extension, '.prelim_n_genes.png')

    # Cannot run filter_cells and filter_genes in backed mode
    adata_mem = adata.to_memory()

    sc.pp.filter_cells(adata_mem, min_genes=3)  # this adds adata.obs.n_genes
    sc.pp.filter_genes(adata_mem, min_cells=300)    # this adds adata.obs.n_cells though we do not use it
    sc.pp.calculate_qc_metrics(adata_mem, inplace=True) # This will get total_counts

    # rename total_counts to n_counts for consistency with the rest of the codebase
    adata_mem.obs['n_counts'] = adata_mem.obs['total_counts']

    sc.pl.violin(adata_mem, ['n_genes', 'n_counts'],
                    jitter=0.4, multi_panel=True, save="_prelim_violin.png")

    sc.pl.scatter(adata_mem, x='n_counts', y='n_genes', save="_prelim_n_genes.png")

    # move files written to tmp
    shutil.move("/tmp/violin_prelim_violin.png", violin_image_path)
    shutil.move("/tmp/scatter_prelim_n_genes.png", scatter_image_path)


def detect_clustering(adata: "AnnData") -> bool:
    """
    Looks first for 'cluster', then 'cell_type'
    """
    cols = adata.obs.columns.tolist()

    for vname in VALID_CLUSTER_COLUMN_NAMES:
        if vname in cols:
            return True

    if has_clustering(adata):
        return True

    return False

def detect_tsne(adata: "AnnData") -> bool:
    """
    Looks for the combination of pairs in VALID_TSNE_PAIRS or existing tSNE in obsm
    """
    cols = adata.obs.columns.tolist()

    for pair in VALID_TSNE_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            return True

    # Maybe it's already in obsm
    if has_tsne(adata):
        return True

    return False

def detect_umap(adata: "AnnData") -> bool:
    """
    Looks for the combination of pairs in VALID_UMAP_PAIRS or existing UMAP in obsm
    """
    cols = adata.obs.columns.tolist()

    for pair in VALID_UMAP_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            return True

    # Maybe it's already in obsm (i.e. spatial uploads)
    if has_umap(adata):
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
        if "X_tsne" in adata.obsm.keys():
            return True
    except Exception:
        pass

    return False

def has_umap(adata: "AnnData") -> bool:
    try:
        if "X_umap" in adata.obsm.keys():
            return True
    except Exception:
        pass

    return False