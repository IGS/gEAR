
import sys
import typing
from pathlib import Path

import colorcet as cc
import numpy as np
import pandas as pd
import spatialdata as sd
from flask import request
from flask_restful import Resource

from .common import create_projection_adata

TWO_LEVELS_UP = 2
ABS_PATH_WWW = Path(__file__).resolve().parents[TWO_LEVELS_UP]  # web-root dir
PANEL_CSV_CACHE_DIR = ABS_PATH_WWW / "cache" / "spatial_panel"

SPATIAL_PATH = ABS_PATH_WWW.joinpath("datasets/spatial")
PROJECTIONS_BASE_DIR = ABS_PATH_WWW.joinpath("projections")

GEAR_ROOT = ABS_PATH_WWW.parent
LIB_PATH = GEAR_ROOT.joinpath("lib")
sys.path.append(str(LIB_PATH))
from gear.spatialhandler import SPATIALTYPE2CLASS

if typing.TYPE_CHECKING:
    from anndata import AnnData
    from gear.spatialhandler import SpatialHandler

def prep_sdata(dataset_id: str) -> "SpatialHandler":
    """
    Prepare and return a SpatialHandler for the given dataset.

    This function locates a .zarr dataset under the module-level SPATIAL_PATH using
    the provided dataset_id, reads it with sd.read_zarr, extracts the platform
    identifier from sdata.tables["table"].uns["platform"], validates that the
    platform is supported by the SPATIALTYPE2CLASS mapping, instantiates the
    appropriate SpatialHandler subclass, attaches the loaded sdata to it and
    returns the handler.

    Parameters
    ----------
    dataset_id : str
        Identifier of the dataset (without the ".zarr" suffix). The function will
        look for a file at SPATIAL_PATH / f"{dataset_id}.zarr".

    Returns
    -------
    SpatialHandler
        An instance of the handler class corresponding to the dataset platform,
        with its `sdata` attribute set to the loaded dataset.

    Raises
    ------
    ValueError
        - If the dataset zarr path does not exist.
        - If the dataset does not contain platform information at
          sdata.tables["table"].uns["platform"].
        - If the platform value is not present in the SPATIALTYPE2CLASS mapping.
        In the case of an unsupported platform, the function will print an error
        and the set of supported types to stderr before raising.

    Notes
    -----
    - This function depends on the module-level names: SPATIAL_PATH, sd,
      SPATIALTYPE2CLASS, and sys.
    - The handler class returned is constructed via SPATIALTYPE2CLASS[platform](),
      and the loaded sdata is assigned to its `sdata` attribute prior to return.
    """
    zarr_path = SPATIAL_PATH / f"{dataset_id}.zarr"
    if not zarr_path.exists():
        raise ValueError(f"Dataset {dataset_id} not found")

    sdata = sd.read_zarr(zarr_path)

    try:
        platform = sdata.tables["table"].uns["platform"]
    except KeyError:
        raise ValueError("No platform information found in the dataset")

    # Ensure the spatial data type is supported
    if platform not in SPATIALTYPE2CLASS.keys():
        print("Invalid or unsupported spatial data type", file=sys.stderr)
        print("Supported types: {0}".format(SPATIALTYPE2CLASS.keys()), file=sys.stderr)
        raise ValueError(
            f"Invalid or unsupported spatial data type: {platform}"
        )

    # Use uploader class to determine correct helper functions
    spatial_obj: "SpatialHandler" = SPATIALTYPE2CLASS[platform]()
    spatial_obj.sdata = sdata
    return spatial_obj

def generate_spatial_image_df(spatial_obj: "SpatialHandler") -> np.ndarray | None:
    """
    Generate a NumPy array representation of the spatial image for a given SpatialHandler.

    This function attempts to extract an image from the provided spatial handler by
    calling its extract_img() method. Some spatial datasets do not include image
    data or may raise errors during extraction/conversion; such errors are caught,
    an error message is printed to sys.stderr, and the function returns None.

    Parameters
    ----------
    spatial_obj : SpatialHandler
        An object that implements an extract_img() method which returns a NumPy array
        representation of the spatial image.

    Returns
    -------
    np.ndarray | None
        A NumPy array containing the image data if extraction succeeds, otherwise
        None when no image is available or an error occurred during extraction.

    Notes
    -----
    - Errors raised by spatial_obj.extract_img() are handled internally; they will
      not be propagated to the caller. Instead, a message describing the failure
      is written to sys.stderr.
    - This function is intended for non-fatal retrieval of optional image data.

    Examples
    --------
    >>> arr = generate_spatial_image_df(spatial_obj)
    >>> if arr is not None:
    ...     # process the NumPy array
    ...     pass
    """

    try:
        return spatial_obj.extract_img()
    except Exception as e:
        # This is not a fatal error. Some spatial datasets do not have images
        print(f"No image found or error converting image to dataframe: {e}", file=sys.stderr)
    return None

def create_gene_df(adata: "AnnData", gene_symbol: str) -> pd.DataFrame:
    """
    Create a pandas DataFrame for a single gene from an AnnData object, including spatial,
    optional UMAP coordinates, and cluster annotations.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing the observations (obs), variables (var),
        and multidimensional annotations (obsm). Must contain a 'gene_symbol'
        column in adata.var and 'spatial1' and 'spatial2' in adata.obs.
    gene_symbol : str
        The gene symbol to select. This will be normalized using the
        normalize_searched_gene(dataset_genes, gene_symbol) helper before selection.

    Returns
    -------
    pandas.DataFrame
        A DataFrame indexed by the observations (cells/spots) for which the gene
        was selected. The DataFrame contains at least the following columns:
          - raw_value: expression values for the requested gene (single-column DataFrame)
          - spatial1, spatial2: spatial coordinates from adata.obs
          - clusters: cluster labels from adata.obs (as a categorical dtype)
          - clusters_cat_codes: integer category codes for clusters (as categorical dtype)
        If UMAP coordinates exist in adata.obsm under the key 'X_umap', the DataFrame
        will also include:
          - UMAP1, UMAP2: first and second UMAP dimensions

    Raises
    ------
    ValueError
        If the requested gene is not found in the dataset after normalization
        or if cluster information is not present in adata.obs.

    Notes
    -----
    - The function uses adata.var['gene_symbol'].unique() to build the searchable
      gene set and resolve the requested gene via normalize_searched_gene.
    - UMAP coordinates are only added if 'X_umap' exists in adata.obsm. If UMAP
      coordinates are missing, they are not computed here (TODO in code).
    - Observations with missing cluster values are dropped from the returned DataFrame.

    Examples
    --------
    >>> df = create_gene_df(adata, "GeneA")
    >>> df.columns
    Index(['raw_value', 'spatial1', 'spatial2', 'UMAP1', 'UMAP2', 'clusters',
           'clusters_cat_codes'], dtype='object')
    """
    adata = adata
    dataset_genes = set(adata.var["gene_symbol"].unique())
    norm_gene_symbol = normalize_searched_gene(dataset_genes, gene_symbol)
    if norm_gene_symbol is None:
        raise ValueError(
            f"Gene '{gene_symbol}' not found in the dataset. Please choose a different gene."
        )
    norm_gene_symbol = norm_gene_symbol
    gene_filter = adata.var.gene_symbol == norm_gene_symbol
    selected = adata[:, gene_filter]
    selected.var.index = pd.Index(["raw_value"])
    dataframe = selected.to_df()

    # Add spatial coords
    dataframe["spatial1"] = selected.obs["spatial1"]
    dataframe["spatial2"] = selected.obs["spatial2"]

    # Add n_genes_by_counts for future filtering if it exists
    dataframe["n_genes_by_counts"] = selected.obs["n_genes_by_counts"]

    # Add UMAP coords if they exist
    if "X_umap" in selected.obsm.keys():
        X, Y = (0, 1)
        dataframe["UMAP1"] = selected.obsm["X_umap"].transpose()[X].tolist()
        dataframe["UMAP2"] = selected.obsm["X_umap"].transpose()[Y].tolist()
    else:
        # TODO: recreate UMAP coords if they do not exist, or if the user wants to recompute them after filtering (which requires a new csv)
        pass

    # Add cluster info
    if "clusters" not in selected.obs:
        raise ValueError("No cluster information found in adata.obs")

    # ! dtypes will not be preserved when writing to CSV and reading back later.
    dataframe["clusters"] = selected.obs["clusters"].astype("category")
    dataframe["clusters_cat_codes"] = dataframe["clusters"].cat.codes.astype("category")

    # Drop any NA clusters
    dataframe = dataframe.dropna(subset=["clusters"])
    return dataframe

def map_colors(dataframe: pd.DataFrame, spatial_img: np.ndarray | None) -> pd.DataFrame:
    # Assuming df is your DataFrame and it has a column "clusters"
    unique_clusters = dataframe["clusters"].unique()
    sorted_clusters = sort_clusters(unique_clusters)

    unique_clusters = sorted_clusters

    if "colors" in dataframe:
        color_map = {
            cluster: dataframe[dataframe["clusters"] == cluster][
                "colors"
            ].to_numpy()[0]
            for cluster in unique_clusters
        }
    else:
        # Some glasbey_bw_colors may not show well on a dark background so use "light" colors if images are not present
        # Prepending b_ to the name will return a list of RGB colors (though glasbey_light seems to already do this)
        swatch_color = (
            cc.b_glasbey_bw if spatial_img is not None else cc.glasbey_light
        )

        color_map = {
            cluster: swatch_color[i % len(swatch_color)]
            for i, cluster in enumerate(unique_clusters)
        }
        # Map the colors to the clusters
        dataframe["colors"] = dataframe["clusters"].map(color_map)
    return dataframe

def sort_clusters(clusters) -> list:
    """
    Sort clusters by number if numerical, otherwise by name.
    """
    try:
        sorted_clusters = sorted(clusters, key=lambda x: int(x))
    except Exception:
        sorted_clusters = sorted(clusters, key=lambda x: str(x))
    return sorted_clusters

def normalize_searched_gene(gene_set, chosen_gene) -> str | None:
    """Convert to case-insensitive version of gene.  Returns None if gene not found in dataset."""
    chosen_gene_lower = chosen_gene.lower()
    for gene in gene_set:
        try:
            if chosen_gene_lower == gene.lower():
                return gene
        except Exception:
            print(gene, file=sys.stderr)
            raise
    return None

class SpatialPanel(Resource):
    """Resource for prepping spatial data to use in Holoviz Panel app.

    Returns
    -------
    File name for the prepared csv file to import in the Panel app.
    """
    def post(self, dataset_id):

        req = request.get_json()
        gene_symbol = req.get('gene_symbol', None)  # gene symbol or projection pattern
        projection_id = req.get('projection_id', None)

        response = {
            "filename": None,
            "success": 0,
            "message": "",
        }

        DATASET_DIR = PANEL_CSV_CACHE_DIR / dataset_id
        DATASET_DIR.mkdir(parents=True, exist_ok=True)

        filename_to_check = f"{gene_symbol}.csv"
        if projection_id:
            filename_to_check = f"{projection_id}_{gene_symbol}.csv"

        # If csv is cached, return it
        csv_path = DATASET_DIR / filename_to_check
        if csv_path.is_file():
            response["filename"] = str(csv_path.name)
            response["message"] = "Using cached file."
            response["success"] = 1
            return response


        try:
            spatial_obj = prep_sdata(dataset_id)
            spatial_img = generate_spatial_image_df(spatial_obj)

            if spatial_img is not None:
                spatial_img_path = DATASET_DIR / "spatial_img.npy"
                np.save(spatial_img_path, spatial_img)

            adata = spatial_obj.sdata.tables["table"]

            # Modify the adata object to use the projection ID if it exists
            if projection_id:
                adata = create_projection_adata(adata, dataset_id, projection_id)

            gene_df = create_gene_df(adata, gene_symbol)
            gene_df = map_colors(gene_df, spatial_img)

            gene_df.to_csv(csv_path, index=False)

            response["success"] = 1
            response["filename"] = str(csv_path.name)

        except Exception as e:
            response["message"] = f"Error preparing data: {e}"
        finally:
            return response

