"""
Module centered around adding a primary analysis to a dataset.  This includes:
- Detecting whether tSNE, UMAP, or clustering analyses have already been performed based
    on the presence of certain columns in the AnnData object
- Applies to single-cell-rnaseq datasets and spatial datasets.
"""

import json
import shutil
from pathlib import Path

import numpy as np
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix, issparse

from .analysis import H5adAdapter, ZarrAdapter

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
    """
    Set up the primary analysis for an uploaded dataset in its staging directory.

    Loads (or creates from template) analysis_pipeline.json, writes preliminary composition
    plots, and records any detected tSNE, UMAP or clustering results, adding them to the
    AnnData object and rewriting the H5AD and JSON files when changes are made.

    Args:
        dataset_id (str): Dataset ID, also used as the primary analysis ID.
        share_id (str): Share ID used to name the uploaded .h5ad or .zarr file.
        staging_dir (Path): Directory containing the uploaded dataset files.
        dataset_format (str): Dataset format; "spatial" loads a .zarr store, anything else an .h5ad.

    Returns:
        bool: True on success.

    Raises:
        PrimaryAnalysisProcessingError: If the uploaded dataset file is not found.
    """

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

    # Ensure Ensembl IDs are not duplicated, which will throw errors downstream
    adata.var_names_make_unique()

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
    """
    Copy the first recognized cluster column in adata.obs into obs['louvain'] as categories.
    """
    cols = adata.obs.columns.tolist()

    for vname in VALID_CLUSTER_COLUMN_NAMES:
        if vname in cols:
            user_defined_cluster_names = adata.obs[vname].astype('category')
            adata.obs['louvain'] = user_defined_cluster_names
            return

def add_tsne_analysis(adata: "AnnData") -> None:
    """
    Copy the first recognized tSNE coordinate column pair in adata.obs into obsm['X_tsne'].
    """
    cols = adata.obs.columns.tolist()

    for pair in VALID_TSNE_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            adata.obsm['X_tsne'] = adata.obs[[pair[0], pair[1]]].values
            return

def add_umap_analysis(adata: "AnnData") -> None:
    """
    Copy the first recognized UMAP coordinate column pair in adata.obs into obsm['X_umap'].
    """
    cols = adata.obs.columns.tolist()

    for pair in VALID_UMAP_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            adata.obsm['X_umap'] = adata.obs[[pair[0], pair[1]]].values
            return

def create_composition_plots(adata: "AnnData", dataset_path: str, is_spatial: bool) -> None:
    """
    Write preliminary QC plots (genes/counts per cell violin and counts-vs-genes scatter).

    Per-cell stats are computed in chunks so backed AnnData objects are never fully
    loaded into memory; cells with fewer than 3 expressed genes are excluded.

    Args:
        adata (AnnData): Dataset to summarize; may be in backed mode.
        dataset_path (str): Dataset path; its .h5ad/.zarr extension is replaced to name the
            ".prelim_violin.png" and ".prelim_n_genes.png" output images.
        is_spatial (bool): Whether the dataset path uses a .zarr extension instead of .h5ad.
    """
    # Create the pathnames for the images
    extension = ".h5ad"
    if is_spatial:
        extension = ".zarr"

    violin_image_path = str(dataset_path).replace(extension, '.prelim_violin.png')
    scatter_image_path = str(dataset_path).replace(extension, '.prelim_n_genes.png')

    # sc.pp.filter_cells/filter_genes refuse backed AnnData, so this used to call
    # adata.to_memory() first -- which for a backed dataset means loading the entire
    # matrix into RAM just for two throwaway QC plots. Compute the same per-cell stats
    # via chunked_X instead, which reads X in row batches and works in backed mode.
    min_genes = 3
    n_genes = np.zeros(adata.n_obs, dtype=np.int64)
    n_counts = np.zeros(adata.n_obs, dtype=np.float64)

    for chunk, start, stop in adata.chunked_X():
        if issparse(chunk):
            n_genes[start:stop] = chunk.getnnz(axis=1)
            n_counts[start:stop] = np.asarray(chunk.sum(axis=1)).ravel()
        else:
            chunk = np.asarray(chunk)
            n_genes[start:stop] = (chunk != 0).sum(axis=1)
            n_counts[start:stop] = chunk.sum(axis=1)

    keep = n_genes >= min_genes
    qc_obs = adata.obs.iloc[keep.nonzero()[0]][[]].copy()
    qc_obs['n_genes'] = n_genes[keep]
    qc_obs['n_counts'] = n_counts[keep]

    # Zero-column shell AnnData: sc.pl.violin/scatter only read .obs for these keys
    qc_adata = AnnData(X=csr_matrix((int(keep.sum()), 0)), obs=qc_obs)

    sc.pl.violin(qc_adata, ['n_genes', 'n_counts'],
                    jitter=0.4, multi_panel=True, save="_prelim_violin.png")

    sc.pl.scatter(qc_adata, x='n_counts', y='n_genes', save="_prelim_n_genes.png")

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
    """
    Return True if adata.obs already has a 'louvain' clustering column.
    """
    cols = adata.obs.columns.tolist()
    # "louvain" is a legacy name when we used the scanpy
    # louvain analysis instead of the modern leiden one.
    if 'louvain' in cols:
        return True
    else:
        return False

def has_tsne(adata: "AnnData") -> bool:
    """
    Return True if adata.obsm already contains 'X_tsne' coordinates.
    """
    try:
        if "X_tsne" in adata.obsm.keys():
            return True
    except Exception:
        pass

    return False

def has_umap(adata: "AnnData") -> bool:
    """
    Return True if adata.obsm already contains 'X_umap' coordinates.
    """
    try:
        if "X_umap" in adata.obsm.keys():
            return True
    except Exception:
        pass

    return False
