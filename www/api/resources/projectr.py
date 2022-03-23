from flask import request
from flask_restful import Resource
import os, json, sys
import geardb
from gear.plotting import get_config
from plotly.utils import PlotlyJSONEncoder

import pandas as pd

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import StrVector

# TODO: Remove and figure out a good pattern directory structure
import plotly.express as px
PATTERN_BASE_DIR = "/var/www/patterns"

def get_analysis(analysis, dataset_id, session_id, analysis_owner_id):
    """Return analysis object based on various factors."""
    # If an analysis is posted we want to read from its h5ad
    if analysis:
        ana = geardb.Analysis(id=analysis['id'], dataset_id=dataset_id,
                                session_id=session_id, user_id=analysis_owner_id)

        try:
            ana.type = analysis['type']
        except:
            user = geardb.get_user_from_session_id(session_id)
            ana.discover_type(current_user_id=user.id)
    else:
        ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5_path = ds.get_file_path()

        # Let's not fail if the file isn't there
        if not os.path.exists(h5_path):
            raise PlotError("No h5 file found for this dataset")
        ana = geardb.Analysis(type='primary', dataset_id=dataset_id)
    return ana

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

def run_projectR_cmd(target_df, loading_df):
    """
    Convert input Pandas dataframes to R matrix.
    Pass the inputs into the projectR function written in R.
    Return Pandas dataframe of the projectR output
    """

    projectR = importr('projectR')

    # Convert from pandas dataframe to R data.frame
    with localconverter(robjects.default_converter + pandas2ri.converter):
       target_r_df = robjects.conversion.py2rpy(target_df)
       loading_r_df = robjects.conversion.py2rpy(loading_df)

    # data.frame to matrix (projectR has no data.frame signature)
    r_matrix = robjects.r["as.matrix"]
    target_r_matrix = r_matrix(target_r_df)
    loading_r_matrix = r_matrix(loading_r_df)
    # Assign Rownames to each matrix
    target_r_matrix.rownames = StrVector(target_df.index)
    loading_r_matrix.rownames = StrVector(loading_df.index)

    # Run project R command.  Get projectionPatterns matrix
    projection_patterns_r_matrix = projectR.projectR(data=target_r_matrix, loadings=loading_r_matrix, full=False)

    # matrix back to data.frame
    r_df = robjects.r["as.data.frame"]
    projection_patterns_r_df = r_df(projection_patterns_r_matrix)

    # Convert from R data.frame to pandas dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
       projection_patterns_df = robjects.conversion.rpy2py(projection_patterns_r_df)

    return projection_patterns_df

class PlotError(Exception):
    """Error based on plotting issues."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)



class ProjectR(Resource):
    """
    ProjectR Container

    """

    def post(self, dataset_id):

        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        analysis = req.get('analysis', None)    # For source "dataset"
        analysis_owner_id = req.get('analysis_owner_id', None)
        #plot_type = req.get('plot_type')
        #gene_symbols = req.get('gene_symbols', [])
        scope = req.get('scope', "pattern")
        input_value = req.get('input_value', None)  # This changes depending on scope
        pattern_value = req.get('pattern_value', 1) # If pattern chosen, use this particular one
        kwargs = req.get("custom_props", {})    # Dictionary of custom properties to use in plot

        # 'dataset_id' is the target dataset to be projected into the pattern space
        try:
            ana = get_analysis(None, dataset_id, session_id, None)
        except PlotError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

        # Using adata with "backed" mode does not work with volcano plot
        adata = ana.get_adata(backed=False)

        adata.obs = order_by_time_point(adata.obs)

        # get a map of all levels for each column
        columns = adata.obs.select_dtypes(include="category").columns.tolist()

        if 'replicate' in columns:
            columns.remove('replicate')

        success = 1
        message = ""

        """
        What has been done
        * A PC was selected as the input pattern

        Steps

        1. Load matrix from DataMTX.tab or Anndata.X.  The observations and genes are largely stripped away
        2. Load ROWmeta DIMRED tab file for a dataset/analysis combo. (pattern file)
        3. Ensure both data matrix and pattern file are transformed so that genes are indexes
        4. Run projectR to get projectionPatterns matrix, which lists relative weights in the matrix dataset.
            * Rows are PCs, Cols are data observations
            * We need to figure out how to re-assign the PCs
        5. Get observation data (COLmeta or Anndata.obs) and split on a conditon or combination of them.
            * One thing is to transpose projectR output,
        6. Extract input pattern row (as index) from projectionPattern
            * NOTE: the current projectionPattern file changes rownames to "x1", "x2", etc., which is different from the ROWmeta DIMRED file
        7. Plot.  This differs depending on plot type, but the example is
          * x - condition groups
          * y - weights per group
          * color - another condition + color mapping if applicable

        Only step 4 needs to be in R, and I guess we could use rpy2 to call that.

        """

        # Ensure target dataset has genes as rows
        target_df = adata.to_df().transpose()

        loading_df = None

        if scope == "repository":
            # Pattern repository
            # Row: Genes
            # Col: Patterns
            loading_df = pd.read_csv("{}/{}".format(PATTERN_BASE_DIR, input_value), sep="\t")

            # Assumes first column is gene info. Standardize on a common index name
            loading_df.rename(columns={ loading_df.columns[0]:"dataRowNames" }, inplace=True)
            loading_df.set_index('dataRowNames', inplace=True)
        elif scope == "dataset":
            # Previous analysis stored in a dataset
            # Row: Genes
            # Col: PCA/tSNE/UMAP 1/2

            # 'dataset_id' is the target dataset to be projected into the pattern space
            try:
                source_ana = get_analysis(analysis, input_value, session_id, analysis_owner_id)
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

            # Using adata with "backed" mode does not work with volcano plot
            source_adata = source_ana.get_adata(backed=False)

            source_adata.obs = order_by_time_point(source_adata.obs)

            # Primary dataset - find tSNE_1 and tSNE_2 in obs and build X_tsne
            if analysis is None or analysis in ["null", "undefined", dataset_id]:
                if "tSNE_1" in source_adata.obs.columns and "tSNE_2" in source_adata.obs.columns:
                    pass
                if "UMAP_1" in source_adata.obs.columns and "UMAP_1" in source_adata.obs.columns:
                    pass
                if "PCA_1" in source_adata.obs.columns and "PCA_2" in source_adata.obs.columns:
                    pass

        elif scope == "gene_cart":
            # Weighted genes present
            # Row; Genes
            # Col: Weights (single col)
            pass
        else:
            raise Exception("Invalid scope was called.")

        # NOTE: This will not work if there are no common genes (i.e. mouse patterns with human dataset)

        projection_patterns_df = run_projectR_cmd(target_df, loading_df).transpose()

        # Append observation categories to list of patterns
        for col in adata.obs.columns:
            projection_patterns_df[col] = adata.obs[col]

        # TODO: Hard-coded - generalize
        # NOTE: Testing with "8353c7d4-ab13-e5db-3675-f03c23a68f7d"
        fig = px.scatter(projection_patterns_df, x="condition", y=pattern_value, color="group")

        # ? What plot types to use?  Single or multi-pattern?
        # ? Pull default display from dataset and sub genes for PCs?

        plot_json = json.dumps(fig, cls=PlotlyJSONEncoder)

        return {
            "success": success
            , "message": message
            , 'plot_json': json.loads(plot_json)
            , "plot_config": get_config()
        }


        return {}