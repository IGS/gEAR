"""
Common functions shared across several resources
"""

import os
import sys
import tempfile
from pathlib import Path

import pandas as pd
from anndata import AnnData
from gear.analysis import (
    Analysis,
    SpatialAnalysis,
    get_analysis,
    normalize_analysis_input,
)
from gear.plotting import PlotError
from shadows import AnnDataShadow
from werkzeug.utils import secure_filename

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP]  # web-root dir
PROJECTIONS_BASE_DIR = abs_path_www.joinpath("projections")

def clip_expression_values(adata: "AnnData", min_clip: float | None=None, max_clip: float | None=None) -> "AnnData":
    """
    Clips the expression values in an AnnData object to specified minimum and/or maximum values.

    Parameters
    ----------
    adata : AnnData
        The AnnData object containing expression data to be clipped.
    min_clip : float or None, optional
        Minimum value to clip the expression data. Values below this will be set to min_clip.
        If None, no minimum clipping is applied.
    max_clip : float or None, optional
        Maximum value to clip the expression data. Values above this will be set to max_clip.
        If None, no maximum clipping is applied.

    Returns
    -------
    AnnData
        The AnnData object with clipped expression values.
    """
    X = adata.to_df()
    X = X.clip(lower=min_clip, upper=max_clip)
    adata.X = X.to_numpy()
    return adata

def get_adata_from_analysis(
    analysis: dict | str | None, dataset_id: str, session_id: str | None, backed: bool = False
) -> AnnData:
    """
    Retrieve an AnnData object associated with a specific analysis.

    Args:
        analysis_id (str | None): The unique identifier of the analysis. Can be None.
        dataset_id (str): The unique identifier of the dataset.
        session_id (str | None): The session identifier. Can be None.

    Returns:
        AnnData: The AnnData object retrieved from the specified analysis.

    Raises:
        Any exceptions raised by `get_analysis` or `ana.get_adata()`.

    """
    analysis_dict = normalize_analysis_input(analysis)
    ana: Analysis = get_analysis(analysis_dict, dataset_id, session_id, is_spatial=False)
    if isinstance(ana, SpatialAnalysis):
        raise ValueError("Analysis is not of type Analysis")
    return ana.get_adata(backed=backed)


def get_adata_shadow_from_analysis(
    analysis: dict | str | None, dataset_id: str, session_id: str | None
) -> AnnDataShadow:
    """
    Retrieve an AnnData object associated with a specific analysis.

    Args:
        analysis_id (dict | str | None): The unique identifier of the analysis. Can be None.
        dataset_id (str): The unique identifier of the dataset.
        session_id (str | None): The session identifier. Can be None.

    Returns:
        AnnDataShadow: An object representing the shadow of the AnnData associated with the analysis.

    Raises:
        Any exceptions raised by get_user_from_session_id, Analysis, or AnnDataShadow constructors.
    """
    analysis_dict = normalize_analysis_input(analysis)
    ana = get_analysis(analysis_dict, dataset_id, session_id)
    return AnnDataShadow(ana.dataset_path)


def get_adata_shadow(
    analysis: dict | str | None,
    dataset_id: str,
    session_id: str | None,
) -> AnnData | AnnDataShadow:
    """
    Retrieve an AnnData or AnnDataShadow object for a given dataset and analysis context.

    Parameters:
        analysis_id (dict | str | None): The analysis identifier, or None if not applicable.
        dataset_id (str): The unique identifier for the dataset.
        session_id (str | None): The session identifier, or None if not applicable.

    Returns:
        AnnData | AnnDataShadow: The loaded AnnData or AnnDataShadow object, depending on the dataset type and context.

    Notes:
        - If `analysis_id` is provided, the AnnDataShadow is loaded from the analysis context; otherwise, from the primary dataset.
    """
    analysis_dict = normalize_analysis_input(analysis)
    return get_adata_shadow_from_analysis(analysis_dict, dataset_id, session_id)

    # see https://github.com/scverse/shadows/issues/4
    # As of 0.1a2, support for legacy AnnData objects is fixed. Leaving this check here for safeguarding.
    # from pandas.api.types import is_integer_dtype
    # if is_integer_dtype(adata.var.index.dtype) :
    #    print("Using AnnData instead of AnnDataShadow because var and obs are not correct for dataset {}".format(dataset_id), file=sys.stderr)
    #    adata = get_adata_from_analysis(analysis_id, dataset_id, session_id)
    # return adata


def get_spatial_adata(
    analysis: dict | str | None,
    dataset_id: str,
    session_id: str | None,
) -> AnnData:
    """
    Retrieve an AnnData object associated with a specific spatial analysis.

    Args:
        analysis_id (dict | str | None): The ID of the analysis to retrieve, or None.
        dataset_id (str): The ID of the dataset to retrieve.
        session_id (str | None): The session ID for user authentication, or None.

    Returns:
        AnnData: The processed spatial AnnData object, with optional image data and metadata.

    Raises:
        ValueError: If the spatial data type is invalid or unsupported.
    """

    # Get spatial-based adata object.
    analysis_dict = normalize_analysis_input(analysis)
    ana = get_analysis(analysis_dict, dataset_id, session_id, is_spatial=True)

    if not isinstance(ana, SpatialAnalysis):
        raise ValueError("Analysis is not spatial")

    adata = ana.get_adata()
    return adata


def create_projection_adata(
    dataset_adata: AnnData, dataset_id: str, projection_id: str
) -> AnnData:
    """
    Creates a new AnnData object representing a projection, using data from a CSV file and metadata from an existing AnnData object.

    Parameters:
        dataset_adata (AnnData): The source AnnData object containing observation (obs), observation metadata (obsm), and unstructured data (uns) to be copied to the new projection.
        dataset_id (str): The identifier for the dataset, used to locate the projection CSV file.
        projection_id (str): The identifier for the projection, used to locate the projection CSV file.

    Returns:
        AnnData: A new AnnData object constructed from the projection CSV file and metadata from the original dataset_adata.

    Raises:
        PlotError: If the projection AnnData object cannot be created from the CSV file.

    Notes:
        - The function reads a CSV file containing the projection data, constructs a new AnnData object with the same obs, obsm, and uns as the original dataset_adata, and sets the gene_symbol in var.
        - The function ensures that file paths are sanitized to prevent path traversal.
        - If dataset_adata is backed, its file is closed after processing to avoid stale file handles.
    """
    # Create AnnData object out of readable CSV file
    projection_id = secure_filename(projection_id)
    dataset_id = secure_filename(dataset_id)

    projection_dir = Path(PROJECTIONS_BASE_DIR).joinpath("by_dataset", dataset_id)
    # Sanitize input to prevent path traversal
    projection_csv_path = projection_dir.joinpath("{}.csv".format(projection_id))
    try:
        import pandas as pd

        # READ CSV to make X and var
        projection_dataframe = pd.read_csv(
            projection_csv_path, sep=",", index_col=0, header=0
        )
        X = projection_dataframe.to_numpy()
        var = pd.DataFrame(index=projection_dataframe.columns)
        obs = dataset_adata.obs
        obsm = dataset_adata.obsm
        uns = dataset_adata.uns

        # This should resolve (https://github.com/IGS/gEAR/issues/951)
        # Create a temporary file
        with tempfile.NamedTemporaryFile(delete=True) as temp_file:
            temp_file_path = temp_file.name  # Associate with the temporary filename to ensure AnnData is read in "backed" mode

            # Create the anndata object in memory mode and write to h5ad (since the dataset_adata is backed)
            projection_adata = AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns) # type: ignore
            # For some reason the gene_symbol is not taken in by the constructor
            projection_adata.var["gene_symbol"] = projection_adata.var_names
            projection_adata.write_h5ad(filename=Path(temp_file_path))

    except Exception as e:
        print(str(e), file=sys.stderr)
        raise PlotError("Could not create projection AnnData object from CSV.")
    finally:
        # Close dataset adata so that we do not have a stale opened object
        if dataset_adata.isbacked:
            dataset_adata.file.close()

    return projection_adata


def order_by_time_point(obs_df: pd.DataFrame) -> pd.DataFrame:
    """
    Reorders the 'time_point' column in a DataFrame based on the 'time_point_order' column, if present.
    If the DataFrame contains a 'time_point_order' column and a 'time_point' column, this function:
        - Sorts the DataFrame by 'time_point_order'.
        - Converts 'time_point' to a categorical type and reorders its categories according to the sorted order.
        - Removes the 'time_point_order' column from the DataFrame.
    Parameters:
        obs_df (pd.DataFrame): Input DataFrame containing at least 'time_point' and optionally 'time_point_order' columns.
    Returns:
        pd.DataFrame: The DataFrame with the 'time_point' column reordered and 'time_point_order' column removed (if present).
    """

    # check if time point order is intially provided in h5ad
    time_point_order = obs_df.get("time_point_order")
    if time_point_order is not None and "time_point" in obs_df.columns:
        sorted_df = obs_df.drop_duplicates().sort_values(by="time_point_order")
        # Safety check. Make sure time point is categorical before
        # calling .cat
        obs_df["time_point"] = pd.Categorical(obs_df["time_point"])
        col = obs_df["time_point"].cat
        obs_df["time_point"] = col.reorder_categories(
            sorted_df.time_point.drop_duplicates(), ordered=True
        )
        obs_df = obs_df.drop(["time_point_order"], axis=1)
    return obs_df
