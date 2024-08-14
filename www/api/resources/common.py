"""
Common functions shared across several resources
"""

import sys
import anndata

from pathlib import Path

from werkzeug.utils import secure_filename

from gear.plotting import PlotError

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP] # web-root dir
PROJECTIONS_BASE_DIR = abs_path_www.joinpath('projections')

def create_projection_adata(dataset_adata, dataset_id, projection_id):
    # Create AnnData object out of readable CSV file
    projection_id = secure_filename(projection_id)
    dataset_id = secure_filename(dataset_id)

    projection_dir = Path(PROJECTIONS_BASE_DIR).joinpath("by_dataset", dataset_id)
    # Sanitize input to prevent path traversal
    projection_adata_path = projection_dir.joinpath("{}.h5ad".format(projection_id))
    projection_csv_path = projection_dir.joinpath("{}.csv".format(projection_id))
    try:
        import pandas as pd
        # READ CSV to make X and var
        df = pd.read_csv(projection_csv_path, sep=',', index_col=0, header=0)
        X = df.to_numpy()
        var = pd.DataFrame(index=df.columns)
        obs = dataset_adata.obs
        # Create the anndata object and write to h5ad
        # Associate with a filename to ensure AnnData is read in "backed" mode
        projection_adata = anndata.AnnData(X=X, obs=obs, var=var, obsm=dataset_adata.obsm, filemode='r')
    except Exception as e:
        print(str(e), file=sys.stderr)
        raise PlotError("Could not create projection AnnData object from CSV.")
    # Close dataset adata so that we do not have a stale opened object
    if dataset_adata.isbacked:
        dataset_adata.file.close()

    # For some reason the gene_symbol is not taken in by the constructor
    projection_adata.var["gene_symbol"] = projection_adata.var_names

    # Associate with a filename to ensure AnnData is read in "backed" mode
    # This creates the h5ad file if it does not exist
    # TODO: If too many processes read from this file, it can throw a BlockingIOError. Eventually we should
    #       handle this by creating a copy of the file for each process, like a tempfile.
    projection_adata.filename = projection_adata_path

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