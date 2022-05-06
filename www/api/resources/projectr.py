from flask import request
from flask_restful import Resource
from pathlib import Path
import os, json, sys
import geardb
import gear.rfuncs as rfx

import pandas as pd

# TODO: Remove and figure out a good pattern directory structure
TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP] # web-root dir
PATTERN_BASE_DIR = os.path.abspath(os.path.join(abs_path_www, 'patterns'))
WEIGHTED_CARTS_BASE_DIR = os.path.abspath(os.path.join(abs_path_www, 'carts'))
PROJECTIONS_BASE_DIR = os.path.abspath(os.path.join(abs_path_www, 'projections'))

VALID_SCOPES = ["projection_pattern", "genecart"]

def build_projectR_output_file(dataset_id, input_value, pca_ext, mkdir=True):
    """Builds the output file path for the projectR analysis."""
    projection_subdir = os.path.join(PROJECTIONS_BASE_DIR, dataset_id)
    if mkdir:
        Path(projection_subdir).mkdir(parents=True, exist_ok=True)  # Create dir if it doesn't exist
    return "{}/{}{}.csv".format(projection_subdir, dataset_id, input_value, pca_ext)



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
            raise ProjectRError("No h5 file found for this dataset")
        ana = geardb.Analysis(type='primary', dataset_id=dataset_id)
    return ana

class ProjectRError(Exception):
    """Error based on issues that would manifest in the projectR call."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)

class ProjectROutputFile(Resource):
    """
    Get or create path to output file. Will mkdir if the directory does not currently exist.
    """

    def post(self, dataset_id):
        req = request.get_json()
        scope = req.get('scope', "pattern")
        input_value = req.get('input_value', None)  # This changes depending on scope
        is_pca = req.get('is_pca', False)

        projection_csv = None

        # Extension to append to filename if PCA analysis is requested
        pca_ext = "_pca" if is_pca else ""

        if scope == "projection_pattern":
            pattern_title = input_value.replace('.tab', '').replace('ROWmeta_DIMRED_', '')
            projection_csv = build_projectR_output_file(dataset_id, pattern_title, pca_ext, mkdir=False)
        elif scope == "genecart":
            input_file = "cart.{}.tab".format(input_value)
            projection_csv = build_projectR_output_file(dataset_id, input_file, pca_ext, mkdir=False)
        else:
            raise Exception("Invalid scope was called.")

        return {
            "file_found": 1 if os.path.exists(projection_csv) else 0,
            "csv_file": projection_csv
        }


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
        scope = req.get('scope', "projection_pattern")
        input_value = req.get('input_value', None)  # This changes depending on scope
        is_pca = req.get('is_pca', False)
        kwargs = req.get("custom_props", {})    # Dictionary of custom properties to use in plot

        # 'dataset_id' is the target dataset to be projected into the pattern space
        try:
            ana = get_analysis(None, dataset_id, session_id, None)
        except ProjectRError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

        # Using adata with "backed" mode does not work with volcano plot
        adata = ana.get_adata(backed=False)

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
        projection_csv = None

        # Extension to append to filename if PCA analysis is requested
        pca_ext = "_pca" if is_pca else ""


        if scope == "projection_pattern":
            # Pattern repository
            # Row: Genes
            # Col: Patterns
            loading_df = pd.read_csv("{}/{}".format(PATTERN_BASE_DIR, input_value), sep="\t")
            pattern_title = input_value.replace('.tab', '').replace('ROWmeta_DIMRED_', '')
            projection_csv = build_projectR_output_file(dataset_id, pattern_title, pca_ext)

        elif scope == "dataset":
            # Previous analysis stored in a dataset
            # Row: Genes
            # Col: PCA/tSNE/UMAP 1/2
            pass
        elif scope == "genecart":
            # Weighted genes present... behaves exactly like the "repository" case
            # Row; Genes
            # Col: Weights
            # The gene cart "share_id" is passed to the CGI script as the file_name, which we can get the actual cart name from
            input_file = "cart.{}.tab".format(input_value)
            loading_df = pd.read_csv("{}/{}".format(WEIGHTED_CARTS_BASE_DIR, input_file), sep="\t")
            projection_csv = build_projectR_output_file(dataset_id, input_file, pca_ext)

        else:
            raise Exception("Invalid scope was called.")

        # Assumes first column is gene info. Standardize on a common index name
        loading_df.rename(columns={ loading_df.columns[0]:"dataRowNames" }, inplace=True)
        loading_df.set_index('dataRowNames', inplace=True)

        # NOTE: This will not work if there are no common genes (i.e. mouse patterns with human dataset)

        # If projectR has already been run, we can just load the csv file.  Otherwise, let it rip!
        if os.path.isfile(projection_csv):
            print("INFO: Found exisitng projection_csv file, loading it.")
        else:

            # Perform overlap to see if there are overlaps between genes from both dataframes
            if target_df.index.intersection(loading_df.index).empty:
                message = "No common genes between the target dataset and the pattern file."
                return {"success": -1, "message": message}

            projection_patterns_df = rfx.run_projectR_cmd(target_df, loading_df, is_pca).transpose()
            projection_patterns_df.to_csv(projection_csv)

        return {
            "success": success
            , "message": message
            , "csv_file": os.path.basename(projection_csv)
        }

