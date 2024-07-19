from flask import request
from flask_restful import Resource, reqparse
from pathlib import Path
import json, sys, fcntl
from io import StringIO
import pandas as pd
import asyncio, aiohttp
import gc

from os import getpid
from time import sleep
from more_itertools import sliced

import geardb
from gear.orthology import get_ortholog_file, map_dataframe_genes

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")
PROJECTIONS_BASE_DIR = abs_path_www.joinpath("projections")
PROJECTIONS_JSON_BASENAME = "projections.json"

ANNOTATION_TYPE = "ensembl" # NOTE: This will change in the future to be varied.

# limit of asynchronous tasks that can happen at a time
# I am setting this slightly under the "MaxKeepAliveRequests" in apache.conf
SEMAPHORE_LIMIT = 50

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
parser.add_argument('algorithm', help="An algorithm needs to be provided", type=str,  required=False)
parser.add_argument('analysis', type=str, required=False)   # not used at the moment

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

def calculate_chunk_size(num_genes, num_samples):
    """
    Calculate number of chunks to divide all samples into.
    """
    TOTAL_DATA_LIMIT = 1e6
    total_data = num_genes * num_samples
    total_data_chunks = total_data / TOTAL_DATA_LIMIT
    return int(num_samples / total_data_chunks) # take floor.  returned value * num_genes < total_data_limit

def concat_fetch_results_to_dataframe(res_jsons):
    # Concatenate the dataframes back together again
    res_dfs = [pd.read_json(StringIO(res_json), orient="split", dtype="float32") for res_json in res_jsons]
    projection_patterns_df = pd.concat(res_dfs)
    return projection_patterns_df

def create_new_uuid(dataset_id, genecart_id, algorithm):
    import hashlib, uuid
    uuid_str = "{}-{}-{}".format(dataset_id, genecart_id, algorithm)
    md5 = hashlib.md5()
    md5.update(uuid_str.encode("utf-8"))
    return uuid.UUID(md5.hexdigest())

def create_unweighted_loading_df(genecart):
    # Now convert into a GeneCollection to get the Ensembl IDs (which will be the unique identifiers)
    gene_collection = geardb.GeneCollection()
    gene_collection.get_by_gene_symbol(gene_symbol=" ".join(genecart.genes), exact=True, organism_id=genecart.organism_id)

    #TODO: Need to filter by organism

    loading_data = ({"dataRowNames": gene.ensembl_id, "gene_sym":gene.gene_symbol, "unweighted":1} for gene in gene_collection.genes)
    return pd.DataFrame(loading_data)

def create_weighted_loading_df(genecart_id):
    file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format("cart." + genecart_id))
    try:
        return pd.read_csv(file_path, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError("Could not find pattern file {}".format(file_path))

def chunk_dataframe(df, chunk_size):
    # Chunk dataset by samples/cells (cols). Is a generator function
    # Help from: https://stackoverflow.com/questions/51674751/using-requests-library-to-make-asynchronous-requests-with-python-3-7
    index_slices = sliced(range(len(df.columns)), chunk_size)
    for idx, index_slice in enumerate(index_slices):
        yield df.iloc[:,index_slice]

def limited_as_completed(coros, limit):
    """A version of asyncio.as_completed that takes a generator of coroutines instead of a list."""
    # Source: https://www.artificialworlds.net/blog/2017/05/31/python-3-large-numbers-of-tasks-with-limited-concurrency/
    # Uses far less memory than the list version (asyncio.as_completed)
    from itertools import islice

    futures = [
        asyncio.ensure_future(c)
        for c in islice(coros, 0, limit)
        ]

    async def first_to_finish():
        while True:
            await asyncio.sleep(0)
            for f in futures:
                if f.done():
                    futures.remove(f)
                    try:
                        newf = next(coros)
                        futures.append(
                            asyncio.ensure_future(newf))
                    except StopIteration as e:
                        pass
                    return f.result()

    while futures:
        yield first_to_finish()

async def fetch_all(target_df, loading_df, algorithm, genecart_id, dataset_id, chunk_size):
    """Create coroutine tasks out of all chunked projection cloud run service POST requests."""

    async with aiohttp.ClientSession() as client:
        # Create coroutines to be executed.
        loadings_json = loading_df.to_json(orient="split")
        coros = (fetch_one(client, {
                "target": chunk_df.to_json(orient="split")
                , "loadings": loadings_json
                , "algorithm": algorithm
                , "genecart_id":genecart_id # This helps in identifying which combinations are going through
                , "dataset_id":dataset_id
                }) for chunk_df in chunk_dataframe(target_df, chunk_size))

        # This loop processes results as they come in.
        return [await res for res in limited_as_completed(coros, SEMAPHORE_LIMIT)]

async def fetch_one(client, payload):
    """
    makes an non-authorized POST request to the specified HTTP endpoint
    """

    # Cloud Run uses your service's hostname as the `audience` value
    # audience = 'https://my-cloud-run-service.run.app/'
    # For Cloud Run, `endpoint` is the URL (hostname + path) receiving the request
    # endpoint = 'https://my-cloud-run-service.run.app/my/awesome/url'

    audience=this.servercfg['projectR_service']['hostname']
    endpoint="{}/".format(audience)
    headers = {"content_type": "application/json"}

    # https://docs.aiohttp.org/en/stable/client_reference.html
    # (semaphore) https://stackoverflow.com/questions/40836800/python-asyncio-semaphore-in-async-await-function
    async with client.post(url=endpoint, json=payload, headers=headers, raise_for_status=True) as response:

        #TODO: Figure out how to get the response status code and payload to print to stderr
        if response.status != 200:
            # print payload to stderr
            print("ERROR: POST request to {} failed with status code {}".format(endpoint, response.status), file=sys.stderr)
            print("ERROR: Payload was: {}".format(payload), file=sys.stderr)

        return await response.json()

def projectr_callback(dataset_id, genecart_id, projection_id, session_id, scope, algorithm, fh):
    success = 1
    message = ""

    if not fh:
        fh=sys.stderr

    if scope == "unweighted-list" and algorithm in ["nmf", "fixednmf"]:
        return {
            'success': -1
            , 'message': "Unweighted gene lists cannot be used with NMF algorithms."
        }

    """
    Steps

    1. Load AnnData object and transform into a dataframe of genes x observations
    2. Load pattern file (gene cart) of genes x patterns. For unweighted gene carts, all weights are 1
    3. Run projectR to get projectionPatterns matrix, which lists relative weights in the matrix dataset.
        * Output - Rows are patterns, Cols are data observations (before transposing)

    Only step 3 needs to be in R, and we use rpy2 to call that.
    """

    # Get unique identifier of first gene from target dataset
    #first_dataset_gene = target_df.index[0]

    # Unweighted carts get a "1" weight for each gene
    genecart = geardb.get_gene_cart_by_share_id(genecart_id)
    if not genecart:
        return {
            'success': -1
            , 'message': "Could not find gene cart in database"
        }

    genecart.get_genes();

    if not len(genecart.genes):
        return {
            'success': -1
            , 'message': "No genes found within this gene cart"
        }

    # Row: Genes
    # Col: Pattern weights
    try:
        loading_df = create_unweighted_loading_df(genecart) if scope == "unweighted-list" else create_weighted_loading_df(genecart_id)
    except Exception as e:
        print(str(e), file=fh)
        import traceback
        traceback.print_exc()
        return {
            'success': -1
            , 'message': str(e)
        }

    # Assumes first column is unique identifiers. Standardize on a common index name
    loading_df.rename(columns={ loading_df.columns[0]:"dataRowNames" }, inplace=True)

    # Drop the gene symbol column
    loading_df = loading_df.drop(loading_df.columns[1], axis=1)

    loading_df.set_index('dataRowNames', inplace=True)

    # If cross-species, remap the genecart genes to the orthologous genes for the dataset's organism
    try:
        # Get the organism ID associated with the dataset.
        try:
            ds = geardb.get_dataset_by_id(dataset_id)
        except:
            raise Exception("Dataset was not found in the database. Please contact a gEAR admin.")

        if not genecart.organism_id == ds.organism_id:
            ortholog_file = get_ortholog_file(genecart.organism_id, ds.organism_id, ANNOTATION_TYPE)
            loading_df = map_dataframe_genes(loading_df, ortholog_file)
    except Exception as e:
        print(str(e), file=fh)
        import traceback
        traceback.print_exc()
        return {"success": -1, "message": str(e)}

    # Drop duplicate unique identifiers. This may happen if two unweighted gene cart genes point to the same Ensembl ID in the db
    loading_df = loading_df[~loading_df.index.duplicated(keep='first')]

    # NOTE Currently no analyses are supported yet.
    try:
        ana = geardb.get_analysis(None, dataset_id, session_id)
    except Exception as e:
        print(str(e), file=fh)
        return {
            'success': -1
            , 'message': str(e)
        }

    # Using adata with "backed" mode does not work with volcano plot
    adata = ana.get_adata(backed=True)

    # If dataset genes have duplicated index names, we need to rename them to avoid errors
    # in collecting rownames in projectR (which gives invalid output)
    # This means these duplicated genes will not be in the intersection of the dataset and pattern genes
    adata.var_names_make_unique()

    num_target_genes = adata.shape[1]
    num_loading_genes = loading_df.shape[0]

    # Perform overlap to see if there are overlaps between genes from both dataframes
    index_intersection = adata.var.index.intersection(loading_df.index)
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

    # Reduce the size of both dataframes before POSTing
    adata = adata[:,index_intersection]
    loading_df = loading_df.loc[index_intersection]

    # Ensure target dataset has genes as rows
    # ! The to_df() seems to be a memory_chokepoint, doubling the memory used by the adata object
    target_df = adata.to_df().transpose()

    # Close dataset adata so that we do not have a stale opened object
    adata.file.close()

    dataset_projection_csv = build_projection_csv_path(dataset_id, projection_id, "dataset")

    # Create lock file if it does not exist
    lockfile = str(dataset_projection_csv) + ".lock"
    if Path(lockfile).exists():
        print("INFO: Found lockfile for another current projectR run of {}.  Going to wait for that run to finish and steal its output.".format(projection_id), file=fh)
        try:
            # Test to see if the exclusive lock has expired
            lock_fh = create_lock_file(lockfile)
            print("INFO: Lock for {} seems to be stale".format(projection_id), file=fh)
            remove_lock_file(lock_fh, lockfile)
        except:
            print("INFO: Lock for {} seems to be valid.".format(projection_id), file=fh)
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

    # Chunk size needs to adjusted by how many genes are present, so that the payload always stays under the body size limit
    chunk_size = calculate_chunk_size(len(target_df.index), len(target_df.columns))
    if len(target_df.columns) < chunk_size:
        chunk_size = len(target_df.columns)

    print("TARGET: {}\nGENECART: {}\nTARGET DF (genes,samples): {}\nSAMPLES PER CHUNK: {}".format(dataset_id, genecart_id, target_df.shape, chunk_size), file=fh)

    # shuffle target dataframe rows.  Needed to balance out the chunks in the NMF algorithms
    target_df = target_df.sample(frac=1)

    if algorithm == "fixednmf":
        # Normalize the target_df by the minimum sum of each expression column
        columns_sums = target_df.sum(axis=0)
        normalized_columns_sums = columns_sums / columns_sums.min()
        target_df = target_df.div(normalized_columns_sums, axis=1)

    if this.servercfg['projectR_service']['cloud_run_enabled'].startswith("1"):
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        try:
            results = loop.run_until_complete(fetch_all(target_df, loading_df, algorithm, genecart_id, dataset_id, chunk_size))
            print("INFO: All fetch tasks have completed", file=fh)
        except Exception as e:
            print(str(e), file=fh)
            # Raises as soon as one "gather" task has an exception
            remove_lock_file(lock_fh, lockfile)
            return {
                'success': -1
                , 'message': "Something went wrong with the projection-creating step."
                , "num_common_genes": intersection_size
                , "num_genecart_genes": num_loading_genes
                , "num_dataset_genes": num_target_genes
            }
        finally:
            # Wait 250 ms for the underlying SSL connections to close
            loop.run_until_complete(asyncio.sleep(0.250))
            loop.stop() # prevent "Task was destroyed but it is pending!" messages
            loop.close()

        print("INFO: Concatenating results to dataframe", file=fh)
        projection_patterns_df = concat_fetch_results_to_dataframe(results)

        del results
        gc.collect()    # trying to clear memory

        if len(projection_patterns_df.index) != len(adata.obs.index):
            message = "Not all chunked sample rows were returned by projectR.  Cannot proceed."
            print(message, file=sys.stderr)
            return {
                'success': -1
                , 'message': message
                , "num_common_genes": intersection_size
                , "num_genecart_genes": num_loading_genes
                , "num_dataset_genes": num_target_genes
            }

        # There is a good chance the samples are now out of order, which will break
        # the copying of the dataset observation metadata when this output is converted
        # to an AnnData object. So reorder back to dataset sample order.
        projection_patterns_df = projection_patterns_df.reindex(adata.obs.index.tolist())
    else:
        # If not using the cloud run service, do this on the server
        abs_path_gear = Path(__file__).resolve().parents[3]
        service_path = abs_path_gear.joinpath('services')
        # abs_path_lib is a Path object so we need to convert to string
        sys.path.insert(0, str(service_path))

        # https://github.com/IGS/gEAR/issues/442#issuecomment-1317239909
        # Basically this is a stopgap until projectR has an option to remove
        # centering around zero for PCA loadings.  Chunking the data breaks
        # the output due to the centering around zero step.
        try:
            if algorithm == "pca":
                from projectr.main import do_pca_projection
                projection_patterns_df = do_pca_projection(target_df,loading_df)
            elif algorithm == "binary":
                from projectr.main import do_binary_projection
                projection_patterns_df = do_binary_projection(target_df,loading_df)
            elif algorithm == "2silca":
                pass
            elif algorithm in ["nmf", "fixednmf"]:
                from projectr.rfuncs import run_projectR_cmd
                projection_patterns_df = run_projectR_cmd(target_df, loading_df, algorithm).transpose()
            else:
                raise ValueError("Algorithm {} is not supported".format(algorithm))
        except Exception as e:
            print(str(e), file=sys.stderr)
            return {
                'success': -1
                , 'message': "Something went wrong with the projection-creating step."
                , "num_common_genes": intersection_size
                , "num_genecart_genes": num_loading_genes
                , "num_dataset_genes": num_target_genes
            }

    # Have had cases where the column names are x1, x2, x3, etc. so load in the original pattern names
    projection_patterns_df = projection_patterns_df.set_axis(loading_df.columns, axis="columns")

    print("INFO: Writing projection patterns to {}".format(dataset_projection_csv), file=fh)
    projection_patterns_df.to_csv(dataset_projection_csv)

    del projection_patterns_df
    gc.collect()    # trying to clear memory

    dataset_projection_json_file = build_projection_json_path(dataset_id, "dataset")

    # Add new configuration to the list for this dictionary key
    with open(dataset_projection_json_file) as projection_fh:
            dataset_projections_dict = json.load(projection_fh)
    dataset_projections_dict.setdefault(genecart_id, []).append({
        "uuid": projection_id
        , "algorithm": algorithm
        , "num_common_genes": intersection_size
        , "num_genecart_genes": num_loading_genes
        , "num_dataset_genes": num_target_genes
    })
    write_to_json(dataset_projections_dict, dataset_projection_json_file)

    # Do the same for the genecart version
    genecart_projection_csv = build_projection_csv_path(genecart_id, projection_id, "genecart")
    genecart_projection_json_file = build_projection_json_path(genecart_id, "genecart")

    # Symlink dataset_projection_csv to genecart_projection_csv (this syntax feels like it's in reverse but it is correct)
    # NOTE: In the Docker instance, symlink reflects path to the mounted volume, not the local path
    try:
        genecart_projection_csv.symlink_to(dataset_projection_csv)
    except FileExistsError:
        print("Symlink already exists for {}".format(dataset_projection_csv), file=fh)

    with open(genecart_projection_json_file) as projection_fh:
        genecart_projections_dict = json.load(projection_fh)
    genecart_projections_dict.setdefault(dataset_id, []).append({
        "uuid": projection_id
        , "algorithm": algorithm
        , "num_common_genes": intersection_size
        , "num_genecart_genes": num_loading_genes
        , "num_dataset_genes": num_target_genes
    })
    write_to_json(genecart_projections_dict, genecart_projection_json_file)

    # Remove file lock
    print("INFO: Removing lock file for {}".format(projection_id), file=fh)
    remove_lock_file(lock_fh, lockfile)

    return {
        "success": success
        , "message": message
        , "projection_id": projection_id
        , "num_common_genes": intersection_size
        , "num_genecart_genes": num_loading_genes
        , "num_dataset_genes": num_target_genes
    }

class ProjectROutputFile(Resource):
    """
    Get or create path to output file. Will mkdir if the directory does not currently exist.
    """

    def post(self, dataset_id):
        args = parser.parse_args()
        genecart_id = args['genecart_id']
        algorithm = args['algorithm']

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

        with open(dataset_projection_json_file) as projection_fh:
            projections_dict = json.load(projection_fh)

        # If the pattern was not projected onto this dataset, initialize a list of configs
        # "cart.{projection_id}" added for backwards compatability
        if not (genecart_id in projections_dict or "cart.{}".format(genecart_id) in projections_dict):
            # NOTE: tried to write JSON with empty list but it seems that empty keys are skipped over.

            return {
                "projection_id": None
            }

        # If legacy version exists, copy to current format.
        if not genecart_id in projections_dict or "cart.{}".format(genecart_id) in projections_dict:
            print("Copying legacy cart.{} to {} in the projection json file.".format(genecart_id, genecart_id), file=fh)
            projections_dict[genecart_id] == projections_dict["cart.{}".format(genecart_id)]
            projections_dict.pop("cart.{}".format(genecart_id), None)
            # move legacy genecart projection stuff to new version
            print("Moving legacy cart.{} contents to {} in the projection genecart directory.".format(genecart_id, genecart_id), file=fh)
            old_genecart_projection_json_file = build_projection_json_path("cart.{}".format(genecart_id), "genecart")
            old_genecart_projection_json_file.parent.rename(dataset_projection_json_file.parent)

        for config in projections_dict[genecart_id]:
            # Handle legacy algorithm
            if "is_pca" in config:
                config["algorithm"] = "pca" if config["is_pca"] == 1 else "nmf"

            if algorithm == config["algorithm"]:
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
    Formats inputs to prep for projectR API call running on Google Cloud Run
    """

    def post(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        args = run_projectr_parser.parse_args()

        genecart_id = args['genecart_id']
        algorithm = args['algorithm']
        projection_id = args['projection_id']
        scope = args['scope']

        projection_id = projection_id or create_new_uuid(dataset_id, genecart_id, algorithm)
        dataset_projection_csv = build_projection_csv_path(dataset_id, projection_id, "dataset")
        dataset_projection_json_file = build_projection_json_path(dataset_id, "dataset")

        # If projectR has already been run, we can just load the csv file.  Otherwise, let it rip!
        if Path(dataset_projection_csv).is_file():
            print("INFO: Found exisitng dataset_projection_csv file {}, loading it.".format(dataset_projection_csv))

            # Projection already exists, so we can just return info we want to return in a message
            with open(dataset_projection_json_file) as projection_fh:
                projections_dict = json.load(projection_fh)
            common_genes = None
            genecart_genes = None
            dataset_genes = None
            for config in projections_dict[genecart_id]:
                # Handle legacy algorithm
                if "is_pca" in config:
                    config["algorithm"] = "pca" if config["is_pca"] == 1 else "nmf"

                if algorithm == config["algorithm"]:
                    common_genes = config.get('num_common_genes', None)
                    genecart_genes = config.get('num_genecart_genes', -1)
                    dataset_genes = config.get('num_dataset_genes', -1)
                    break

            if common_genes:
                message = "Found {} common genes between the target dataset ({} genes) and the pattern file ({} genes).".format(common_genes, dataset_genes, genecart_genes)

            return {
                "success": 1
                , "message": message
                , "projection_id": projection_id
                , "num_common_genes": common_genes
                , "num_genecart_genes": genecart_genes
                , "num_dataset_genes": dataset_genes
            }


        # Create a messaging queue if necessary. Make it persistent across the lifetime of the Flask server.
        # Channels will be spawned during each task.
        if this.servercfg['projectR_service']['queue_enabled'].startswith("1"):
            import gearqueue
            host = this.servercfg['projectR_service']['queue_host']
            try:
                # Connect as a blocking RabbitMQ publisher
                connection = gearqueue.Connection(host=host, publisher_or_consumer="publisher")
            except Exception as e:
                return {
                    'success': -1
                , 'message': str(e)
                }
            # Connect as a blocking RabbitMQ publisher
            with connection:
                connection.open_channel()
                task_finished = False
                response = {}
                def _on_response(channel, method_frame, properties, body):
                    nonlocal task_finished
                    nonlocal response
                    task_finished = True
                    response = json.loads(body)
                    print("[x] - Received response for dataset {} and genecart {}".format(payload["dataset_id"], payload["genecart_id"]), file=sys.stderr)

                # Create a "reply-to" consumer
                # see https://pika.readthedocs.io/en/stable/examples/direct_reply_to.html?highlight=reply_to#direct-reply-to-example
                try:
                    connection.replyto_consume(
                        on_message_callback=_on_response
                    )
                except Exception as e:
                    return {
                        'success': -1
                    , 'message': str(e)
                    }

                # Create the publisher
                payload = dict()
                payload["dataset_id"] = dataset_id
                payload["session_id"] = session_id
                payload["projection_id"] = projection_id    # the "arg" was assigned and potentially modified
                payload["genecart_id"] = genecart_id
                payload["scope"] = scope
                payload["algorithm"] = algorithm

                try:
                    connection.publish(
                        queue_name="projectr"
                        , reply_to="amq.rabbitmq.reply-to"
                        , message=payload   # method dumps JSON
                    )
                    print("[x] Requesting for dataset {} and genecart {}".format(dataset_id, genecart_id), file=sys.stderr)
                except Exception as e:
                    return {
                        'success': -1
                    , 'message': str(e)
                    }
                # Wait for callback to finish, then return the response
                while not task_finished:
                    pass
                print("[x] sending payload response back to client for dataset {} and genecart {}".format(dataset_id, genecart_id), file=sys.stderr)
                return response
        else:
            return projectr_callback(dataset_id, genecart_id, projection_id, session_id, scope, algorithm, sys.stderr)
