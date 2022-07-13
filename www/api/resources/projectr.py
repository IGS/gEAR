from flask import request
from flask_restful import Resource, reqparse
from pathlib import Path
import json, uuid
import geardb
import gear.rfuncs as rfx
from gear.rfuncs import RError

import pandas as pd

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")
PROJECTIONS_BASE_DIR = abs_path_www.joinpath("projections")
PROJECTIONS_JSON_BASENAME = "projections.json"

"""
projections json format - one in each "projections/by_dataset/<dataset_id> subdirectory

{
  <cart.abc123>: [
    configuration options dict * N configs
    ],
  <cart.def456: [
    configuration_options dict * N configs
    ]
}

Also one in each "projections/by_genecart/<genecart_id> subdirectory.
Total number of projections in whole by_genecart directory = total number in by_dataset directory.

{
    <dataset123>: [
        configuration options dict * N configs
    ],
    <dataset456: [
        configuration options dict * N configs
    ],
}

"""

parser = reqparse.RequestParser(bundle_errors=True)
parser.add_argument('genecart_id', help='Weighted (pattern) genecart id required', type=str, required=True)
parser.add_argument('is_pca', help="'is_pca' needs to be a boolean", type=bool,  required=False)
parser.add_argument('analysis', type=str, required=False)   # not used at the moment
parser.add_argument('analysis_owner_id', type=str, required=False)  # Not used at the moment

run_projectr_parser = parser.copy()
run_projectr_parser.add_argument('projection_id', type=str, required=False)

def build_projection_csv_path(dir_id, file_id, scope):
    """Build the path to the csv file for a given projection. Returns a Path object."""
    return Path(PROJECTIONS_BASE_DIR).joinpath("by_{}".format(scope), dir_id, "{}.csv".format(file_id))

def build_projection_json_path(dir_id, scope):
    """Build the path to the projections json for a given dataset or genecart directory. Returns a Path object."""
    return Path(PROJECTIONS_BASE_DIR).joinpath("by_{}".format(scope), dir_id, PROJECTIONS_JSON_BASENAME)

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
            raise FileNotFoundError("No h5 file found for this dataset")
        ana = geardb.Analysis(type='primary', dataset_id=dataset_id)
    return ana



class ProjectROutputFile(Resource):
    """
    Get or create path to output file. Will mkdir if the directory does not currently exist.
    """

    def post(self, dataset_id):
        args = parser.parse_args()
        genecart_id = args['genecart_id']
        is_pca = args['is_pca']

        # Create the directory if it doesn't exist
        # NOTE: The "mkdir" and "touch" commands are not atomic. Do not fail if directory or file was created by another process.
        dataset_projection_json_file = build_projection_json_path(dataset_id, "dataset")
        if not dataset_projection_json_file.parent.is_dir():
            dataset_projection_json_file.parent.mkdir(parents=True, exist_ok=True)

        genecart_projection_json_file = build_projection_json_path(genecart_id, "genecart")
        if not genecart_projection_json_file.parent.is_dir():
            genecart_projection_json_file.parent.mkdir(parents=True, exist_ok=True)

        # Create the genecart projection json file if it doesn't exist
        # We will use the dataset json file as a check if a projection already exists though
        if not genecart_projection_json_file.is_file():
            genecart_projection_json_file.touch(exist_ok=True)
            write_to_json({}, genecart_projection_json_file) # create empty json file

        # If the dataset json file exists, we can read it and get the output file path
        if not dataset_projection_json_file.is_file():
            dataset_projection_json_file.touch(exist_ok=True)
            write_to_json({}, dataset_projection_json_file) # create empty json file
            return {
                "projection_id": None
            }

        projections_dict = json.load(open(dataset_projection_json_file))

        # If the pattern was not projected onto this dataset, initialize a list of configs
        if not genecart_id in projections_dict:
            # NOTE: tried to write JSON with empty list but it seems that empty keys are skipped over.
            return {
                "projection_id": None
            }

        for config in projections_dict[genecart_id]:
            if int(is_pca) == config['is_pca']:
                projection_id = config['uuid']
                dataset_projection_csv = build_projection_csv_path(dataset_id, projection_id, "dataset")
                return {
                    "projection_id": projection_id  if Path(dataset_projection_csv).is_file() else None
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

        success = 1
        message = ""

        genecart_id = args['genecart_id']
        is_pca = args['is_pca']
        output_id = args['projection_id']

        projection_id = output_id if output_id else uuid.uuid4()
        dataset_projection_csv = build_projection_csv_path(dataset_id, projection_id, "dataset")
        dataset_projection_json_file = build_projection_json_path(dataset_id, "dataset")

        genecart_projection_csv = build_projection_csv_path(genecart_id, projection_id, "genecart")
        genecart_projection_json_file = build_projection_json_path(genecart_id, "genecart")

        # If projectR has already been run, we can just load the csv file.  Otherwise, let it rip!
        if Path(dataset_projection_csv).is_file():
            print("INFO: Found exisitng dataset_projection_csv file {}, loading it.".format(dataset_projection_csv))

            # Projection already exists, so we can just return info we want to return in a message
            projections_dict = json.load(open(dataset_projection_json_file))
            common_genes = None
            genecart_genes = None
            dataset_genes = None
            for config in projections_dict[genecart_id]:
                if int(is_pca) == config['is_pca']:
                    common_genes = config.get('num_common_genes', None)
                    genecart_genes = config.get('num_genecart_genes', -1)
                    dataset_genes = config.get('num_dataset_genes', -1)
                    break

            if common_genes:
                message = "Found {} common genes between the target dataset and the pattern file.".format(common_genes)

            return {
                "success": success
                , "message": message
                , "projection_id": projection_id
                , "num_common_genes": common_genes
                , "num_genecart_genes": genecart_genes
                , "num_dataset_genes": dataset_genes
            }

        """
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

        Only step 4 needs to be in R, and we use rpy2 to call that.
        """

        # NOTE Currently no analyses are supported yet.
        try:
            ana = get_analysis(None, dataset_id, session_id, None)
        except Exception as e:
            return {
                'success': -1
                , 'message': str(e)
            }

        # Using adata with "backed" mode does not work with volcano plot
        adata = ana.get_adata(backed=False)

        # Ensure target dataset has genes as rows
        target_df = adata.to_df().transpose()
        loading_df = None

        # Row: Genes
        # Col: Pattern weights
        # TODO: prioritize reading from h5ad file.
        file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format(genecart_id))
        try:
            loading_df = pd.read_csv(file_path, sep="\t")
        except FileNotFoundError:
            return {
                'success': -1
                , 'message': "Could not find pattern file {}".format(file_path)
            }

        # Store gene symbol series before dropping later
        #gene_syms_series = loading_df[1]

        # Assumes first column is unique identifiers. Standardize on a common index name
        loading_df.rename(columns={ loading_df.columns[0]:"dataRowNames" }, inplace=True)

        # Drop the gene symbol column
        loading_df = loading_df.drop(loading_df.columns[1], axis=1)
        loading_df.set_index('dataRowNames', inplace=True)

        num_target_genes = target_df.shape[0]
        num_loading_genes = loading_df.shape[0]

        # Perform overlap to see if there are overlaps between genes from both dataframes
        index_intersection = target_df.index.intersection(loading_df.index)
        intersection_size = index_intersection.size

        # If no overlap on the genecart unique identifiers, try to overlap on the dropped gene symbols column
        #if index_intersection.empty:
        #    series_intersection = target_df.index.intersection(gene_syms_series)
        #    intersection_size = series_intersection.size

        # If gene symbols overlap, replace loading_df index with gene_syms_series
        #if intersection_size:
        #    pass


        if intersection_size == 0:
            message = "No common genes between the target dataset ({} genes) and the pattern file ({} genes).".format(num_target_genes, num_loading_genes)
            return {
                "success": -1
                , "message": message
                , "num_common_genes": intersection_size
                , "num_genecart_genes": num_loading_genes
                , "num_dataset_genes": num_target_genes
            }

        message = "Found {} common genes between the target dataset ({} genes) and the pattern file ({} genes).".format(intersection_size, num_target_genes, num_loading_genes)

        """
        # NOTE: This is not needed for now, but may happen later
        # Perform a mapping of the dataset genes to the loading file genes
        # 1. Get dataset db entry using dataset_id, including the organism ID
        # 2. For every gene in the loading file, determine if they overlap with the dataset adata.var.index
        # 3. If they do not, then perform a db search for each gene in the loading file to get an Ensembl ID for the common organism
        # 4. Re-map the genes in the loading file to the Ensembl IDs
        ds = geardb.get_dataset_by_id(dataset_id)
        try:
            organism_id = ds.organism_id
        except:
            message = "Dataset was not found in the database. Please contact a gEAR admin."
            return {"success": -1, "message": message}
        """

        try:
            projection_patterns_df = rfx.run_projectR_cmd(target_df, loading_df, is_pca).transpose()
        except RError as re:
            return {
                'success': -1
                , 'message': str(re)
                , "num_common_genes": intersection_size
                , "num_genecart_genes": num_loading_genes
                , "num_dataset_genes": num_target_genes
            }

        # Have had cases where the column names are x1, x2, x3, etc. so load in the original pattern names
        projection_patterns_df.set_axis(loading_df.columns, axis="columns", inplace=True)
        projection_patterns_df.to_csv(dataset_projection_csv)

        # Symlink dataset_projection_csv to genecart_projection_csv (this syntax feels like it's in reverse)
        # NOTE: In the Docker instance, symlink reflects path to the mounted volume, not the local path
        genecart_projection_csv.symlink_to(dataset_projection_csv)

        # Add new configuration to the list for this dictionary key
        dataset_projections_dict = json.load(open(dataset_projection_json_file))
        dataset_projections_dict.setdefault(genecart_id, []).append({
            "uuid": projection_id
            , "is_pca": int(is_pca)
            , "num_common_genes": intersection_size
            , "num_genecart_genes": num_loading_genes
            , "num_dataset_genes": num_target_genes
        })
        write_to_json(dataset_projections_dict, dataset_projection_json_file)

        # Do the same for the genecart version
        genecart_projections_dict = json.load(open(genecart_projection_json_file))
        genecart_projections_dict.setdefault(dataset_id, []).append({
            "uuid": projection_id
            , "is_pca": int(is_pca)
            , "num_common_genes": intersection_size
            , "num_genecart_genes": num_loading_genes
            , "num_dataset_genes": num_target_genes
        })
        write_to_json(genecart_projections_dict, genecart_projection_json_file)

        return {
            "success": success
            , "message": message
            , "projection_id": projection_id
            , "num_common_genes": intersection_size
            , "num_genecart_genes": num_loading_genes
            , "num_dataset_genes": num_target_genes
        }


