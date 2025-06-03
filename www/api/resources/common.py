"""
Common functions shared across several resources
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import anndata
import pandas as pd
from gear.plotting import PlotError
from geardb import Analysis, SpatialAnalysis, get_user_from_session_id
from pandas.api.types import is_integer_dtype
from shadows import AnnDataShadow
from werkzeug.utils import secure_filename

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP] # web-root dir
PROJECTIONS_BASE_DIR = abs_path_www.joinpath('projections')

def get_adata_from_analysis(analysis_id, dataset_id, session_id):
    user = get_user_from_session_id(session_id)
    user_id = None
    if user:
        user_id = user.id
    ana = Analysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user_id)
    ana.discover_type()
    return ana.get_adata()

def get_adata_shadow_from_analysis(analysis_id, dataset_id, session_id):
    user = get_user_from_session_id(session_id)
    user_id = None
    if user:
        user_id = user.id
    ana = Analysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user_id)
    ana.discover_type()
    return AnnDataShadow(ana.dataset_path())

def get_adata_shadow_from_primary(h5_path):
    if not os.path.exists(h5_path):
        raise FileNotFoundError("No h5 file found for this dataset")
    return AnnDataShadow(h5_path)

def get_adata_shadow(analysis_id, dataset_id, session_id, dataset_path, include_images=None):
    if dataset_path.endswith(".zarr"):
        # It's spatial data... probably can't use AnnDataShadow
        return get_spatial_adata(analysis_id, dataset_id, session_id, include_images=include_images)

    if analysis_id:
        adata = get_adata_shadow_from_analysis(analysis_id, dataset_id, session_id)
    else:
        adata = get_adata_shadow_from_primary(dataset_path)

    # see https://github.com/scverse/shadows/issues/4
    # As of 0.1a2, support for legacy AnnData objects is fixed. Leaving this check here for safeguarding.
    if is_integer_dtype(adata.var.index.dtype) :
        print("Using AnnData instead of AnnDataShadow because var and obs are not correct for dataset {}".format(dataset_id), file=sys.stderr)
        adata = get_adata_from_analysis(analysis_id, dataset_id, session_id)
    return adata

def get_spatial_adata(analysis_id: str | None, dataset_id: str, session_id: str, include_images: bool|None = None) -> anndata.AnnData:
    # Get spatial-based adata object.
    user = get_user_from_session_id(session_id)
    user_id = None
    if user:
        user_id = user.id
    ana = SpatialAnalysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user_id)
    ana.discover_type()
    sdata = ana.get_sdata()
    platform = ana.determine_platform(sdata)

    from gear import spatial_handler

    # Ensure the spatial data type is supported
    if not platform or platform not in spatial_handler.SPATIALTYPE2CLASS.keys():
        raise ValueError("Invalid or unsupported spatial data type {0}".format(platform))

    spatial_obj = spatial_handler.SPATIALTYPE2CLASS[platform]()
    spatial_obj.sdata = sdata

    # Filter by bounding box (mostly for images)
    spatial_obj._filter_sdata_by_coords()

    if include_images is None:
        include_images = spatial_obj.has_images

    # Create AnnData object
    # Do not include images in the adata object (to make it lighter)
    spatial_obj._convert_sdata_to_adata(include_images=include_images)
    adata = spatial_obj.adata
    # Extra metadata to help with determining if images are available and where they are
    adata.uns["has_images"] = spatial_obj.has_images
    if spatial_obj.has_images:
        adata.uns["img_name"] = spatial_obj.img_name
    return adata

def create_projection_adata(dataset_adata, dataset_id, projection_id):
    # Create AnnData object out of readable CSV file
    projection_id = secure_filename(projection_id)
    dataset_id = secure_filename(dataset_id)

    projection_dir = Path(PROJECTIONS_BASE_DIR).joinpath("by_dataset", dataset_id)
    # Sanitize input to prevent path traversal
    projection_csv_path = projection_dir.joinpath("{}.csv".format(projection_id))
    try:
        import pandas as pd
        # READ CSV to make X and var
        projection_dataframe = pd.read_csv(projection_csv_path, sep=',', index_col=0, header=0)
        X = projection_dataframe.to_numpy()
        var = pd.DataFrame(index=projection_dataframe.columns)
        obs = dataset_adata.obs
        obsm = dataset_adata.obsm
        uns = dataset_adata.uns

        # This should resolve (https://github.com/IGS/gEAR/issues/951)
        # Create a temporary file
        with tempfile.NamedTemporaryFile(delete=True) as temp_file:
            temp_file_path = temp_file.name      # Associate with the temporary filename to ensure AnnData is read in "backed" mode

            # Create the anndata object in memory mode and write to h5ad (since the dataset_adata is backed)
            projection_adata = anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns)
            # For some reason the gene_symbol is not taken in by the constructor
            projection_adata.var["gene_symbol"] = projection_adata.var_names
            projection_adata.write_h5ad(filename=Path(temp_file_path))
            print(projection_adata.isbacked, file=sys.stderr)

    except Exception as e:
        print(str(e), file=sys.stderr)
        raise PlotError("Could not create projection AnnData object from CSV.")
    finally:
        # Close dataset adata so that we do not have a stale opened object
        if dataset_adata.isbacked:
            dataset_adata.file.close()

    return projection_adata

def order_by_time_point(obs_df):
    """Order observations by time point column if it exists."""
    # check if time point order is intially provided in h5ad
    time_point_order = obs_df.get('time_point_order')
    if (time_point_order is not None and 'time_point' in obs_df.columns):
        sorted_df = obs_df.drop_duplicates().sort_values(by='time_point_order')
        # Safety check. Make sure time point is categorical before
        # calling .cat
        obs_df['time_point'] = pd.Categorical(obs_df['time_point'])
        col = obs_df['time_point'].cat
        obs_df['time_point'] = col.reorder_categories(
            sorted_df.time_point.drop_duplicates(), ordered=True)
        obs_df = obs_df.drop(['time_point_order'], axis=1)
    return obs_df