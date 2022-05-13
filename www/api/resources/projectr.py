from flask import request
from flask_restful import Resource, reqparse
from pathlib import Path
import json, sys, uuid
import geardb
import gear.rfuncs as rfx

import pandas as pd

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")
PROJECTIONS_BASE_DIR = abs_path_www.joinpath("projections")
PROJECTIONS_JSON_BASENAME = "projections.json"

"""
projections json format - one in each "projections/<dataset_id> subdirectory

{
  <pattern_source>: [
    configuration options dict * N configs
    ],
  <pattern_source 2>: [
    configuration_options dict * N configs
    ]
}

"""

parser = reqparse.RequestParser(bundle_errors=True)
parser.add_argument('source_id', help='Pattern source name required', type=str, required=True)
parser.add_argument('is_pca', help="'is_pca' needs to be a boolean", type=bool,  required=False)
parser.add_argument('analysis', type=str, required=False)   # not used at the moment
parser.add_argument('analysis_owner_id', type=str, required=False)  # Not used at the moment

run_projectr_parser = parser.copy()
run_projectr_parser.add_argument('projection_id', type=str, required=False)

def build_projection_csv_path(dataset_id, projection_id):
    return Path(PROJECTIONS_BASE_DIR).joinpath(dataset_id, "{}.csv".format(projection_id))

def build_projection_json_path(dataset_id):
    return Path(PROJECTIONS_BASE_DIR).joinpath(dataset_id, PROJECTIONS_JSON_BASENAME)

def write_to_json(projections_dict, projection_json_file):
    with open(projection_json_file, 'w') as f:
        json.dump(projections_dict, f, ensure_ascii=False, indent=4)


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
        if not Path(h5_path).is_file():
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
        args = parser.parse_args()

        if not Path(PROJECTIONS_BASE_DIR).joinpath(dataset_id).is_dir():
            Path(PROJECTIONS_BASE_DIR).joinpath(dataset_id).mkdir(parents=True)

        projection_json_file = build_projection_json_path(dataset_id)

        # If the json file exists, we can read it and get the output file path
        if not Path(projection_json_file).is_file():
            Path(projection_json_file).touch()
            write_to_json({}, projection_json_file) # create empty json file
            return {
                "projection_id": None
            }

        source_id = args['source_id']
        is_pca = args['is_pca']

        projections_dict = json.load(open(projection_json_file))

        # If the pattern was not projected onto this dataset, initialize a list of configs
        if not source_id in projections_dict:
            # NOTE: tried to write JSON with empty list but it seems that empty keys are skipped over.
            return {
                "projection_id": None
            }

        for config in projections_dict[source_id]:
            if int(is_pca) == config['is_pca']:
                projection_id = config['uuid']
                projection_csv = build_projection_csv_path(dataset_id, projection_id)
                return {
                    "projection_id": projection_id  if Path(projection_csv).is_file() else None
                }

        # If we get here, we didn't find a file for this config
        return {
            "projection_id": None
        }


class ProjectR(Resource):
    """
    ProjectR Container

    """

    def post(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        args = run_projectr_parser.parse_args()

        # Currently no analyses are supported yet.
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
        2. Load ROWmeta DIMRED tab file for a dataset/analysis combo. (pattern file). For weighted gene carts, the format is the same.
        3. Ensure both data matrix and pattern file are transformed so that genes are indexes
        4. Run projectR to get projectionPatterns matrix, which lists relative weights in the matrix dataset.
            * Rows are patterns, Cols are data observations
            * We need to figure out how to re-assign the patterns
        5. Get observation data (COLmeta or Anndata.obs) and split on a conditon or combination of them.
            * One thing is to transpose projectR output,
        6. Extract input pattern row (as index) from projectionPattern

        Only step 4 needs to be in R, and I guess we could use rpy2 to call that.

        """

        source_id = args['source_id']
        is_pca = args['is_pca']
        output_id = args['projection_id']

        # Ensure target dataset has genes as rows
        target_df = adata.to_df().transpose()
        loading_df = None
        projection_id = output_id if output_id else uuid.uuid4()

        projection_csv = build_projection_csv_path(dataset_id, projection_id)

        # Row: Genes
        # Col: Pattern weights
        file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format(source_id))
        loading_df = pd.read_csv(file_path, sep="\t")

        # Assumes first column is gene info. Standardize on a common index name
        loading_df.rename(columns={ loading_df.columns[0]:"dataRowNames" }, inplace=True)
        loading_df.set_index('dataRowNames', inplace=True)

        # NOTE: This will not work if there are no common genes (i.e. mouse patterns with human dataset)

        # If projectR has already been run, we can just load the csv file.  Otherwise, let it rip!
        if Path(projection_csv).is_file():
            print("INFO: Found exisitng projection_csv file {}, loading it.".format(projection_csv))
        else:

            # Perform overlap to see if there are overlaps between genes from both dataframes
            if target_df.index.intersection(loading_df.index).empty:
                message = "No common genes between the target dataset and the pattern file."
                return {"success": -1, "message": message}

            projection_patterns_df = rfx.run_projectR_cmd(target_df, loading_df, is_pca).transpose()
            # Have had cases where the column names are x1, x2, x3, etc. so load in the original pattern names
            projection_patterns_df.set_axis(loading_df.columns, axis="columns", inplace=True)
            projection_patterns_df.to_csv(projection_csv)

            # Add new configuration to the list for this dictionary key
            projection_json_file = build_projection_json_path(dataset_id)
            projections_dict = json.load(open(projection_json_file))
            projections_dict.setdefault(source_id, []).append({
                "uuid": projection_id
                 , "is_pca": int(is_pca)
            })
            write_to_json(projections_dict, projection_json_file)


        return {
            "success": success
            , "message": message
            , "projection_id": projection_id
        }

