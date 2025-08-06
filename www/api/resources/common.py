"""
Common functions shared across several resources
"""

import os
import sys
import tempfile
from pathlib import Path

import pandas as pd
from anndata import AnnData
from gear.plotting import PlotError
from geardb import Analysis, SpatialAnalysis, get_user_from_session_id
from shadows import AnnDataShadow
from werkzeug.utils import secure_filename

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP]  # web-root dir
PROJECTIONS_BASE_DIR = abs_path_www.joinpath("projections")


def get_adata_from_analysis(
    analysis_id: str | None, dataset_id: str, session_id: str | None
) -> AnnData:
    """
    Retrieve an AnnData object associated with a given analysis.

    Args:
        analysis_id (str | None): The unique identifier for the analysis. Can be None.
        dataset_id (str): The unique identifier for the dataset.
        session_id (str | None): The session identifier for the user. Can be None.

    Returns:
        AnnData: The AnnData object corresponding to the specified analysis.

    Raises:
        Any exceptions raised by `get_user_from_session_id`, `Analysis`, or `ana.get_adata()`.

    Notes:
        - If `session_id` is provided and valid, the user ID is associated with the analysis.
        - The function initializes an `Analysis` object, determines its type, and retrieves the associated AnnData.
    """
    user = get_user_from_session_id(session_id)
    user_id = None
    if user:
        user_id = user.id
    ana = Analysis(
        id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user_id
    )
    if ana.type is None:
        ana.discover_type()
    return ana.get_adata()


def get_adata_shadow_from_analysis(
    analysis_id: str | None, dataset_id: str, session_id: str | None
) -> AnnDataShadow:
    """
    Retrieve an AnnDataShadow object based on analysis, dataset, and session identifiers.

    This function constructs an Analysis object using the provided analysis ID, dataset ID, and session ID.
    It attempts to resolve the user from the session ID, if available, and associates the user ID with the analysis.
    After discovering the analysis type, it returns an AnnDataShadow object corresponding to the dataset path.

    Args:
        analysis_id (str | None): The unique identifier for the analysis, or None.
        dataset_id (str): The unique identifier for the dataset.
        session_id (str | None): The session identifier, or None.

    Returns:
        AnnDataShadow: An object representing the shadow of the AnnData associated with the analysis.

    Raises:
        Any exceptions raised by get_user_from_session_id, Analysis, or AnnDataShadow constructors.
    """
    user = get_user_from_session_id(session_id)
    user_id = None
    if user:
        user_id = user.id
    ana = Analysis(
        id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user_id
    )
    if ana.type is None:
        ana.discover_type()
    return AnnDataShadow(ana.dataset_path())


def get_adata_shadow_from_primary(h5_path: str) -> "AnnDataShadow":
    """
    Creates and returns an AnnDataShadow object from the specified HDF5 file path.

    Args:
        h5_path (str): The file path to the HDF5 (.h5) file.

    Returns:
        AnnDataShadow: An instance of AnnDataShadow initialized with the given file path.

    Raises:
        FileNotFoundError: If the specified HDF5 file does not exist at the given path.
    """
    if not os.path.exists(h5_path):
        raise FileNotFoundError("No h5 file found for this dataset")
    return AnnDataShadow(h5_path)


def get_adata_shadow(
    analysis_id: str | None,
    dataset_id: str,
    session_id: str | None,
    dataset_path: str,
    include_images: bool | None = None,
) -> AnnData | AnnDataShadow:
    """
    Retrieve an AnnData or AnnDataShadow object for a given dataset and analysis context.

    Parameters:
        analysis_id (str | None): The analysis identifier, or None if not applicable.
        dataset_id (str): The unique identifier for the dataset.
        session_id (str | None): The session identifier, or None if not applicable.
        dataset_path (str): The file path to the dataset.
        include_images (bool | None, optional): Whether to include images for spatial data. Defaults to None.

    Returns:
        AnnData | AnnDataShadow: The loaded AnnData or AnnDataShadow object, depending on the dataset type and context.

    Notes:
        - If the dataset path ends with ".zarr", spatial data is assumed and a spatial AnnData object is returned.
        - If `analysis_id` is provided, the AnnDataShadow is loaded from the analysis context; otherwise, from the primary dataset.
        - For legacy AnnData objects with integer indices in `.var`, a standard AnnData object is returned for compatibility.
    """
    if dataset_path.endswith(".zarr"):
        # It's spatial data... probably can't use AnnDataShadow
        return get_spatial_adata(
            analysis_id, dataset_id, session_id, include_images=include_images
        )

    if analysis_id:
        return get_adata_shadow_from_analysis(analysis_id, dataset_id, session_id)
    else:
        return get_adata_shadow_from_primary(dataset_path)

    # see https://github.com/scverse/shadows/issues/4
    # As of 0.1a2, support for legacy AnnData objects is fixed. Leaving this check here for safeguarding.
    # from pandas.api.types import is_integer_dtype
    # if is_integer_dtype(adata.var.index.dtype) :
    #    print("Using AnnData instead of AnnDataShadow because var and obs are not correct for dataset {}".format(dataset_id), file=sys.stderr)
    #    adata = get_adata_from_analysis(analysis_id, dataset_id, session_id)
    # return adata


def get_spatial_adata(
    analysis_id: str | None,
    dataset_id: str,
    session_id: str | None,
    include_images: bool | None = None,
) -> AnnData:
    """
    Retrieve a spatial AnnData object for a given analysis and dataset.

    This function initializes a SpatialAnalysis object using the provided analysis ID, dataset ID, and session ID.
    It determines the spatial data platform, validates its support, and processes the spatial data accordingly.
    Optionally, it includes image data in the resulting AnnData object.

    Args:
        analysis_id (str | None): The ID of the analysis to retrieve, or None.
        dataset_id (str): The ID of the dataset to retrieve.
        session_id (str | None): The session ID for user authentication, or None.
        include_images (bool | None, optional): Whether to include images in the AnnData object.
            If None, includes images if available.

    Returns:
        AnnData: The processed spatial AnnData object, with optional image data and metadata.

    Raises:
        ValueError: If the spatial data type is invalid or unsupported.
    """
    # Get spatial-based adata object.
    user = get_user_from_session_id(session_id)
    user_id = None
    if user:
        user_id = user.id
    ana = SpatialAnalysis(
        id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user_id
    )
    if ana.type is None:
        ana.discover_type()
    sdata = ana.get_sdata()
    platform = ana.determine_platform(sdata)

    from gear import spatialhandler

    # Ensure the spatial data type is supported
    if not platform or platform not in spatialhandler.SPATIALTYPE2CLASS.keys():
        raise ValueError(
            "Invalid or unsupported spatial data type {0}".format(platform)
        )

    spatial_obj = spatialhandler.SPATIALTYPE2CLASS[platform]()
    spatial_obj.sdata = sdata

    # Filter by bounding box (mostly for images)
    spatial_obj.filter_sdata_by_coords()

    if include_images is None:
        include_images = spatial_obj.has_images

    # Create AnnData object
    # Do not include images in the adata object (to make it lighter)
    spatial_obj.convert_sdata_to_adata(include_images=include_images)
    adata = spatial_obj.adata
    # Extra metadata to help with determining if images are available and where they are
    adata.uns["has_images"] = spatial_obj.has_images
    if spatial_obj.has_images:
        adata.uns["img_name"] = spatial_obj.img_name
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
            projection_adata = AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns)
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
