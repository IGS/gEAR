from flask import request
from flask_restful import Resource, reqparse
from pathlib import Path
import json, hashlib, uuid, sys, fcntl
import pandas as pd
import requests

from os import getpid
from time import sleep

import geardb

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")
PROJECTIONS_BASE_DIR = abs_path_www.joinpath("projections")
ORTHOLOG_BASE_DIR = abs_path_www.joinpath("feature_mapping")
PROJECTIONS_JSON_BASENAME = "projections.json"

ANNOTATION_TYPE = "ensembl" # NOTE: This will change in the future to be varied.

"""
projections json format - one in each "projections/by_dataset/<dataset_id> subdirectory

{
  <abc123>: [
    configuration options dict * N configs
    ],
  <def456: [
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
parser.add_argument('scope', type=str, required=False)
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

def create_lock_file(filepath):
    """Create an exclusive lock file for the given filepath.  Return the file descriptor."""
    fd = open(filepath, "w+")
    fd.write("{}\n".format(getpid()))
    fcntl.flock(fd, fcntl.LOCK_EX)
    return fd

def remove_lock_file(fd, filepath):
    """Release the lock file."""
    #fcntl.flock(fd, fcntl.LOCK_UN)
    fd.close()
    Path(filepath).unlink()

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

def remap_df_genes(orig_df: pd.DataFrame, orthomap_file: str):
    """Remap the passed-in Dataframe to have gene indexes from the orthologous mapping file."""
    # Read HDF5 file using Pandas read_hdf
    try:
        orthomap_df = pd.read_hdf(orthomap_file)
    except Exception as e:
        print(str(e), file=sys.stderr)
        raise
    # Index -> gs1 / id2 / gs2
    orthomap_dict = orthomap_df.to_dict()["id2"]
    # NOTE: Not all genes can be mapped. Unmappable genes do not change in the original dataframe.
    return orig_df.rename(index=orthomap_dict)

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
        # "cart.{projection_id}" added for backwards compatability
        if not (genecart_id in projections_dict or "cart.{}".format(genecart_id) in projections_dict):
            # NOTE: tried to write JSON with empty list but it seems that empty keys are skipped over.

            return {
                "projection_id": None
            }

        # If legacy version exists, copy to current format.
        if not genecart_id in projections_dict or "cart.{}".format(genecart_id) in projections_dict:
            print("Copying legacy cart.{} to {} in the projection json file.".format(genecart_id, genecart_id), file=sys.stderr)
            projections_dict[genecart_id] == projections_dict["cart.{}".format(genecart_id)]
            projections_dict.pop("cart.{}".format(genecart_id), None)
            # move legacy genecart projection stuff to new version
            print("Moving legacy cart.{} contents to {} in the projection genecart directory.".format(genecart_id, genecart_id), file=sys.stderr)
            old_genecart_projection_json_file = build_projection_json_path("cart.{}".format(genecart_id), "genecart")
            old_genecart_projection_json_file.parent.rename(dataset_projection_json_file.parent)

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
        scope = args['scope']

        uuid_str = "{}-{}-{}".format(dataset_id, genecart_id, is_pca)
        md5 = hashlib.md5()
        md5.update(uuid_str.encode("utf-8"))
        new_uuid = uuid.UUID(md5.hexdigest())

        projection_id = output_id if output_id else new_uuid
        dataset_projection_csv = build_projection_csv_path(dataset_id, projection_id, "dataset")
        dataset_projection_json_file = build_projection_json_path(dataset_id, "dataset")

        genecart_projection_csv = build_projection_csv_path(genecart_id, projection_id, "genecart")
        genecart_projection_json_file = build_projection_json_path(genecart_id, "genecart")

        # If projectR has already been run, we can just load the csv file.  Otherwise, let it rip!
        if Path(dataset_projection_csv).is_file():
            print("INFO: Found exisitng dataset_projection_csv file {}, loading it.".format(dataset_projection_csv), file=sys.stderr)

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

        1. Load AnnData object and transform into a dataframe of genes x observations
        2. Load pattern file (gene cart) of genes x patterns. For unweighted gene carts, all weights are 1
        3. Run projectR to get projectionPatterns matrix, which lists relative weights in the matrix dataset.
            * Output - Rows are patterns, Cols are data observations

        Only step 3 needs to be in R, and we use rpy2 to call that.
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

        # If dataset genes have duplicated index names, we need to rename them to avoid errors
        # in collecting rownames in projectR (which gives invalid output)
        # This means these duplicated genes will not be in the intersection of the dataset and pattern genes
        adata.var_names_make_unique()

        # Ensure target dataset has genes as rows
        target_df = adata.to_df().transpose()

        # Row: Genes
        # Col: Pattern weights
        loading_df = None

        # Get the organism ID associated with the dataset.
        ds = geardb.get_dataset_by_id(dataset_id)
        try:
            dataset_organism_id = ds.organism_id
        except:
            message = "Dataset was not found in the database. Please contact a gEAR admin."
            return {"success": -1, "message": message}

        # Get unique identifier of first gene from target dataset
        #first_dataset_gene = target_df.index[0]

        # Unweighted carts get a "1" weight for each gene
        gc = geardb.get_gene_cart_by_share_id(genecart_id)
        if not gc:
            return {
                'success': -1
                , 'message': "Could not find gene cart in database"
            }

        if not len(gc.genes):
            return {
                'success': -1
                , 'message': "No genes found within this gene cart"
            }

        genecart_organism_id = gc.organism_id

        if scope == "unweighted-list":
            # Now convert into a GeneCollection to get the Ensembl IDs (which will be the unique identifiers)
            gene_collection = geardb.GeneCollection()
            gene_collection.get_by_gene_symbol(gene_symbol=" ".join(gc.genes), exact=True)

            loading_data = []
            for gene in gene_collection.genes:
                loading_data.append({"dataRowNames": gene.ensembl_id, "gene_sym":gene.gene_symbol, "unweighted":1})
            loading_df = pd.DataFrame(loading_data)
        else:
            # weighted carts reside in tab files
            # TODO: prioritize reading from h5ad file.
            file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format("cart." + genecart_id))
            try:
                loading_df = pd.read_csv(file_path, sep="\t")
            except FileNotFoundError:
                return {
                    'success': -1
                    , 'message': "Could not find pattern file {}".format(file_path)
                }

        # Assumes first column is unique identifiers. Standardize on a common index name
        loading_df.rename(columns={ loading_df.columns[0]:"dataRowNames" }, inplace=True)

        # Drop the gene symbol column
        loading_df = loading_df.drop(loading_df.columns[1], axis=1)

        loading_df.set_index('dataRowNames', inplace=True)
        # Drop duplicate unique identifiers. This may happen if two unweighted gene cart genes point to the same Ensembl ID in the db
        loading_df = loading_df[~loading_df.index.duplicated(keep='first')]

        # Get unique identifier of first gene from loading genecart
        #first_loading_gene = loading_df.index[0]

        # If cross-species, remap the genecart genes to the orthologous genes for the dataset's organism
        if not genecart_organism_id == dataset_organism_id:
            orthomap_file_base = "orthomap.{0}.{2}__{1}.{2}.hdf5".format(genecart_organism_id, dataset_organism_id, ANNOTATION_TYPE)
            orthomap_file = ORTHOLOG_BASE_DIR.joinpath(orthomap_file_base)
            try:
                loading_df = remap_df_genes(loading_df, orthomap_file)
            except:
                message = "Could not remap pattern genes to ortholog equivalent in the dataset"
                return {
                    "success": -1
                    , "message": message
                }

        num_target_genes = target_df.shape[0]
        num_loading_genes = loading_df.shape[0]

        # Perform overlap to see if there are overlaps between genes from both dataframes
        index_intersection = target_df.index.intersection(loading_df.index)
        intersection_size = index_intersection.size

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

        # Create lock file if it does not exist
        lockfile = str(dataset_projection_csv) + ".lock"
        if Path(lockfile).exists():
            print("INFO: Found lockfile for another current projectR run of {}.  Going to wait for that run to finish and steal its output.".format(projection_id), file=sys.stderr)
            try:
                # Test to see if the exclusive lock has expired
                lock_fh = create_lock_file(lockfile)
                print("INFO: Lock for {} seems to be stale".format(projection_id), file=sys.stderr)
                remove_lock_file(lock_fh, lockfile)
            except:
                print("INFO: Lock for {} seems to be valid.".format(projection_id), file=sys.stderr)
                # If lock belongs to a valid run, wait for lock to be removed,
                # then return info that is normally returned after projectR is run
                while True:
                    sleep(1)
                    if not Path(lockfile).exists():
                        return {
                            "success": 2
                            , "message": message
                            , "projection_id": projection_id
                            , "num_common_genes": intersection_size
                            , "num_genecart_genes": num_loading_genes
                            , "num_dataset_genes": num_target_genes
                        }

        try:
            lock_fh = create_lock_file(lockfile)
            # NOTE: Will not trigger if process crashes or is forcibly killed off.
        except IOError:
            # This should ideally never be encountered as the previous code should handle existing locked files
            message = "Could not create lock file for this projectR run."
            return {
                'success': -1
                , 'message': message
                , "num_common_genes": intersection_size
                , "num_genecart_genes": num_loading_genes
                , "num_dataset_genes": num_target_genes
            }


        projectr_payload = {
            "target": json.loads(target_df.to_json())
            , "loading": json.loads(loading_df.to_json())
            , "is_pca": is_pca
        }

        response = requests.post(url="TEST", data=projectr_payload, headers=['content_type': 'application.json'])
        if response.raise_for_status():
            # Remove file lock
            remove_lock_file(lock_fh, lockfile)
            return {
                'success': -1
                , 'message': "Could not run projectR command. Job returned a status code of {}".format(response.status_code)
                , "num_common_genes": intersection_size
                , "num_genecart_genes": num_loading_genes
                , "num_dataset_genes": num_target_genes
            }
        projection_patterns_df = pd.read_json(json.dumps(response.json))

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

        # Remove file lock
        remove_lock_file(lock_fh, lockfile)

        return {
            "success": success
            , "message": message
            , "projection_id": projection_id
            , "num_common_genes": intersection_size
            , "num_genecart_genes": num_loading_genes
            , "num_dataset_genes": num_target_genes
        }
