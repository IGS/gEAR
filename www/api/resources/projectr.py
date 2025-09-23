import asyncio
import fcntl
import functools
import gc
import hashlib
import json
import sys
import traceback
import uuid
from io import StringIO
from os import getpid
from pathlib import Path
from time import sleep
from typing import TextIO

import aiohttp
import geardb
import pandas as pd
import scipy.stats as stats
from aiohttp_retry import ExponentialRetry, RetryClient
from flask import request
from flask_restful import Resource, reqparse
from gear.analysis import get_analysis, SpatialAnalysis
from gear.orthology import get_ortholog_file, map_dataframe_genes
from gear.utils import catch_memory_error
from more_itertools import sliced
from werkzeug.utils import secure_filename

# Have all print statements flush immediately (for debugging)
print = functools.partial(print, flush=True)

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig  # noqa: E402

this.servercfg = ServerConfig().parse()  # type: ignore

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP]  # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")
PROJECTIONS_BASE_DIR = abs_path_www.joinpath("projections")
JOB_STATUS_DIR = Path(PROJECTIONS_BASE_DIR).joinpath("job_status")
CHUNK_OUTPUTS_DIR = Path(PROJECTIONS_BASE_DIR).joinpath("chunk_outputs")
PROJECTIONS_JSON_BASENAME = "projections.json"

ANNOTATION_TYPE = "ensembl"  # NOTE: This will change in the future to be varied.

# limit of asynchronous tasks that can happen at a time
# I am setting this slightly under the "MaxKeepAliveRequests" in apache.conf
CONCURRENT_REQUEST_LIMIT = 75
# timeout for the POST request to projectR service
REQUEST_TIMEOUT = 1200

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

"""

parser = reqparse.RequestParser(bundle_errors=True)
parser.add_argument(
    "genecart_id",
    help="Weighted (pattern) genecart id required",
    type=str,
    required=True,
)
parser.add_argument("scope", type=str, required=False)
parser.add_argument(
    "algorithm", help="An algorithm needs to be provided", type=str, required=False
)
parser.add_argument(
    "zscore",
    help="If true, compute z-score calculation before running projectR (or assume it has been calculated)",
    type=bool,
    default=False,
    required=False,
)
parser.add_argument(
    "full_output",
    help="If true, return the full output of the projection, which includes a p-value matrix",
    type=bool,
    default=True,
    required=False,
)
parser.add_argument("analysis", type=str, required=False)  # not used at the moment

run_projectr_parser = parser.copy()
run_projectr_parser.add_argument("projection_id", type=str, required=False)


def build_projection_csv_path(dir_id: str, file_id: str, scope: str) -> Path:
    """Build the path to the csv file for a given projection. Returns a Path object."""
    if scope == "pval":
        # pval files are extra output for the standard "dataset" projections
        return Path(PROJECTIONS_BASE_DIR).joinpath(
            "by_dataset", dir_id, "{}_pval.csv".format(file_id)
        )

    return Path(PROJECTIONS_BASE_DIR).joinpath(
        "by_{}".format(scope), dir_id, "{}.csv".format(file_id)
    )


def build_projection_json_path(dir_id: str, scope: str) -> Path:
    """Build the path to the projections json for a given dataset or genecart directory. Returns a Path object."""
    return Path(PROJECTIONS_BASE_DIR).joinpath(
        "by_{}".format(scope), dir_id, PROJECTIONS_JSON_BASENAME
    )


def create_lock_file(filepath: str) -> TextIO:
    """Create an exclusive lock file for the given filepath.  Return the file descriptor."""
    fd = open(filepath, "w+")
    fd.write("{}\n".format(getpid()))
    fcntl.flock(fd, fcntl.LOCK_EX)
    return fd


def remove_lock_file(fd: TextIO, filepath: str) -> None:
    """Release the lock file."""
    # fcntl.flock(fd, fcntl.LOCK_UN)
    fd.close()
    try:
        Path(filepath).unlink()
    except FileNotFoundError:
        # This is fine, as the lock file may have been removed by another process
        pass


def write_to_json(projections_dict: dict, projection_json_file: Path) -> None:
    with open(projection_json_file, "w") as f:
        json.dump(projections_dict, f, ensure_ascii=False, indent=4)


def calculate_chunk_size(num_genes: int, num_samples: int) -> int:
    """
    Calculate number of chunks to divide all samples into.
    """
    TOTAL_DATA_LIMIT = 1e6
    total_data = num_genes * num_samples
    total_data_chunks = total_data / TOTAL_DATA_LIMIT
    return int(
        num_samples / total_data_chunks
    )  # take floor.  returned value * num_genes < total_data_limit


def concat_fetch_results_to_dataframe(
    projection_id: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Concatenate the dataframes back together again for both "projection" and "pval" keys.
    "pval" may be an empty dataframe.
    """
    projection_dfs = []
    pval_dfs = []

    # Find all the chunk output files for this projection_id
    # If some chunked jobs failed, we will check after this function is called and re-run those chunks
    for filepath in CHUNK_OUTPUTS_DIR.glob(f"{projection_id}_chunk*.json"):
        with open(filepath, "r") as f:
            res_json = json.load(f)
            # res_json is a dictionary. Each value is a JSON string
            if "projection" in res_json:
                projection_json = res_json["projection"]
                projection_df = pd.read_json(StringIO(projection_json), orient="split")
                projection_dfs.append(projection_df)
            if "pval" in res_json:
                pval_json = res_json["pval"]
                pval_df = pd.read_json(StringIO(pval_json), orient="split")
                pval_dfs.append(pval_df)

    projection_patterns_df = (
        pd.concat(projection_dfs) if projection_dfs else pd.DataFrame()
    )
    pval_patterns_df = pd.concat(pval_dfs) if pval_dfs else pd.DataFrame()

    return projection_patterns_df, pval_patterns_df


def create_new_uuid(*args) -> uuid.UUID:
    """
    Generates a new UUID based on the provided arguments.

    Args:
        *args: Variable length argument list. Each argument will be converted to a string and concatenated with a hyphen.

    Returns:
        uuid.UUID: A UUID object generated from the MD5 hash of the concatenated string of arguments.
    """

    uuid_str = "-".join(map(str, args))
    md5 = hashlib.md5()
    md5.update(uuid_str.encode("utf-8"))
    return uuid.UUID(md5.hexdigest())


def create_unweighted_loading_df(genecart: geardb.GeneCart) -> pd.DataFrame:
    # Now convert into a GeneCollection to get the Ensembl IDs (which will be the unique identifiers)
    gene_collection = geardb.GeneCollection()
    genecart.get_genes()
    gene_collection.get_by_gene_symbol(
        gene_symbol=" ".join(genecart.genes),
        exact=True,
        organism_id=genecart.organism_id,
    )

    # TODO: Need to filter by organism

    loading_data = (
        {"dataRowNames": gene.ensembl_id, "gene_sym": gene.gene_symbol, "unweighted": 1}
        for gene in gene_collection.genes
    )
    return pd.DataFrame(loading_data)


def create_weighted_loading_df(genecart_id: str) -> pd.DataFrame:
    file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format("cart." + genecart_id))
    try:
        return pd.read_csv(file_path, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError("Could not find pattern file {}".format(file_path))


def chunk_dataframe(df: pd.DataFrame, chunk_size: int, fh: TextIO):
    # Chunk dataset by samples/cells (cols). Is a generator function
    # Help from: https://stackoverflow.com/questions/51674751/using-requests-library-to-make-asynchronous-requests-with-python-3-7
    index_slices = sliced(range(len(df.columns)), chunk_size)

    for idx, index_slice in enumerate(index_slices):
        yield idx, df.iloc[:, list(index_slice)]

def init_job_status(projection_id: str) -> dict:
    return {"status": "pending", "result": {"projection_id":projection_id}, "error": None}

def write_result_to_file(result, filename) -> None:
    # Write chunked dataframe results to a file, using the projection ID and the dataframe indexes in the filename
    filepath = CHUNK_OUTPUTS_DIR.joinpath(filename)
    with open(filepath, "w") as f:
        json.dump(result, f, indent=4)  # indent for debugging

async def fetch_all_queue(
    target_df: pd.DataFrame,
    loading_df: pd.DataFrame,
    algorithm: str,
    full_output: bool,
    projection_id: str,
    chunk_size: int,
    fh: TextIO,
    concurrency: int = CONCURRENT_REQUEST_LIMIT,
) -> None:
    """
    Asynchronously processes a DataFrame in chunks using a producer-consumer pattern with concurrency control.

    This function splits the `target_df` DataFrame into chunks, enqueues processing tasks, and uses multiple worker coroutines to process each chunk concurrently. Each chunk is sent as a payload to an external service via HTTP requests, and the results are written to disk. The function ensures that already-processed chunks are skipped and supports retry logic for failed requests.

    Args:
        target_df (pd.DataFrame): The DataFrame containing the target data to be processed in chunks.
        loading_df (pd.DataFrame): The DataFrame containing loading data to be included in each payload.
        algorithm (str): The algorithm identifier to be used in the payload.
        full_output (bool): Whether to request full output from the external service.
        projection_id (str): Unique identifier for the projection, used for output file naming.
        chunk_size (int): The number of rows per chunk.
        fh (TextIO): File handle for logging progress and errors.
        concurrency (int, optional): The number of concurrent worker coroutines. Defaults to CONCURRENT_REQUEST_LIMIT.

    Returns:
        None

    Raises:
        Exception: If any worker task encounters an unhandled exception.

    Side Effects:
        - Writes result files for each processed chunk to the output directory.
        - Logs progress and errors to the provided file handle.
    """

    loadings_json = loading_df.to_json(orient="split")
    queue = asyncio.Queue(maxsize=concurrency * 2)

    #total_chunks = sum(1 for _ in chunk_dataframe(target_df, chunk_size, fh))

    async with aiohttp.ClientSession() as client:
        retry_options = ExponentialRetry(
            attempts=3, start_timeout=0.5, max_timeout=0, statuses={500, 502, 503, 504}
        )
        async with RetryClient(client_session=client, retry_options=retry_options, raise_for_status=True) as retry_client:
            # Producer: puts coroutines into the queue
            async def producer():
                for chunk_idx, chunk_df in chunk_dataframe(target_df, chunk_size, fh):

                    # ? Should I add startcol and endcol indexes as well
                    filename = f"{projection_id}_chunk{chunk_idx}.json"
                    filepath = CHUNK_OUTPUTS_DIR.joinpath(filename)
                    if filepath.is_file():
                        print(f"{projection_id} - Chunk {chunk_idx} already processed, skipping.", flush=True, file=fh)
                        continue

                    payload = {
                        "target": chunk_df.to_json(orient="split"),
                        "loadings": loadings_json,
                        "algorithm": algorithm,
                        "full_output": full_output,
                        "projection_id": projection_id,
                        "chunk_idx": chunk_idx,
                    }
                    await queue.put((retry_client, payload, fh))
                # Signal to workers that production is done (sentinel value). One for each worker
                for _ in range(concurrency):
                    await queue.put(None)

            # Worker: gets tasks from the queue and awaits them
            async def worker(idx: int):
                #print(f"{dataset_id} - Worker {idx} started.", flush=True, file=fh)
                try:
                    while True:
                        try:
                            item = await queue.get()
                            if item is None:
                                #print(f"{dataset_id} - Worker {idx} received sentinel value, exiting.", flush=True, file=fh)
                                break
                            #print(f"{dataset_id} - Worker {idx} processing job. Remaining queue size: {queue.qsize()}", flush=True, file=fh)
                            retry_client, payload, filehandle = item
                            try:
                                result = await fetch_one(retry_client, payload, filehandle)
                                if "chunk_idx" not in payload:
                                    raise KeyError("chunk_idx missing from payload")
                                chunk_idx = payload["chunk_idx"]
                                chunk_filename = f"{projection_id}_chunk{chunk_idx}.json"
                                write_result_to_file(result, chunk_filename)
                            except Exception as e:
                                print(f"{projection_id} - Worker {idx} encountered an error: {e}", flush=True, file=fh)
                                print(traceback.format_exc(), file=fh)
                        finally:
                            queue.task_done()
                except Exception as e:
                    print(f"{projection_id} - Worker {idx} crashed with exception: {e}", flush=True, file=fh)
                    print(traceback.format_exc(), flush=True, file=fh)

            # Start producer and workers
            producer_task = asyncio.create_task(producer())
            worker_tasks = [asyncio.create_task(worker(idx)) for idx in range(concurrency)]

            await producer_task
            await queue.join()
            try:
                for w in worker_tasks:
                    if w.done() and w.exception():
                        print(f"Worker task exception: {w.exception()}", flush=True, file=fh)
                    await w
                print(f"{projection_id} - All worker tasks completed successfully.", flush=True, file=fh)
            except Exception as e:
                print(f"{projection_id} - Error in worker tasks: {e}", flush=True, file=fh)
                raise Exception(f"Error in worker tasks: {e}") from e

async def fetch_one(client: RetryClient, payload: dict, fh: TextIO) -> dict:
    """
    makes an non-authorized POST request to the specified HTTP endpoint
    """

    # Cloud Run uses your service's hostname as the `audience` value
    # audience = 'https://my-cloud-run-service.run.app/'
    # For Cloud Run, `endpoint` is the URL (hostname + path) receiving the request
    # endpoint = 'https://my-cloud-run-service.run.app/my/awesome/url'

    audience = this.servercfg["projectR_service"]["hostname"]
    endpoint = "{}/".format(audience)
    headers = {"content_type": "application/json"}

    # https://docs.aiohttp.org/en/stable/client_reference.html
    # (semaphore) https://stackoverflow.com/questions/40836800/python-asyncio-semaphore-in-async-await-function

    dataset_id = payload.get("dataset_id", "unknown")

    try:
        timeout = aiohttp.ClientTimeout(total=REQUEST_TIMEOUT)
        async with client.post(
            url=endpoint, json=payload, headers=headers, timeout=timeout
        ) as response:
            return await response.json()
    except aiohttp.ClientResponseError as cre:
        print(f"{dataset_id} - ERROR: POST request failed with status code {cre.status}", file=fh)
        print(f"{dataset_id} - ERROR: Response body: {cre.message}", file=fh)
        raise aiohttp.ClientResponseError(
            status=cre.status,
            message=f"POST request failed with status code {cre.status}: {cre.message}",
            headers=cre.headers,
            request_info=cre.request_info,
            history=cre.history,
        ) from cre
    except aiohttp.ClientError as ce:
        print(f"{dataset_id} - ERROR: Client error occurred: {str(ce)}", file=fh)
        raise aiohttp.ClientError(
            f"Client error occurred: {str(ce)}",
        ) from ce
    except asyncio.TimeoutError as te:
        print(f"{dataset_id} - ERROR: POST request timed out", file=fh)
        raise asyncio.TimeoutError(
            f"POST request to {endpoint} timed out after {REQUEST_TIMEOUT} seconds"
        ) from te

def write_projection_status(file, status):
    with open(file, "w") as fh:
        json.dump(status, fh)

@catch_memory_error()
def projectr_callback(
    dataset_id: str,
    genecart_id: str,
    projection_id: str,
    session_id: str,
    scope: str,
    algorithm: str,
    zscore: bool,
    full_output: bool,
    fh: TextIO,
) -> dict:
    success = 1
    message = ""

    if not zscore:
        zscore = False

    if not full_output:
        full_output = False

    if not fh:
        fh = sys.stderr

    status = {"status": "pending", "result": {"projection_id": projection_id}, "error": None}
    JOB_STATUS_FILE = JOB_STATUS_DIR.joinpath(f"job_{projection_id}.json")

    if scope == "unweighted-list" and algorithm in ["nmf", "fixednmf"]:
        status["status"] = "failed"
        status["error"] = "Unweighted gene lists cannot be used with NMF algorithms."
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    """
    Steps

    1. Load AnnData object and transform into a dataframe of genes x observations
    2. Load pattern file (gene cart) of genes x patterns. For unweighted gene carts, all weights are 1
    3. Run projectR to get projectionPatterns matrix, which lists relative weights in the matrix dataset.
        * Output - Rows are patterns, Cols are data observations (before transposing)

    Only step 3 needs to be in R, and we use rpy2 to call that.
    """

    # Unweighted carts get a "1" weight for each gene
    genecart = geardb.get_gene_cart_by_share_id(genecart_id)
    if not genecart:
        status["status"] = "failed"
        status["error"] = "Could not find gene list in database."
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    genecart.get_gene_counts()

    if not genecart.num_genes:
        status["status"] = "failed"
        status["error"] = "No genes found within this gene list."
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    # Row: Genes
    # Col: Pattern weights
    try:
        loading_df = (
            create_unweighted_loading_df(genecart)
            if scope == "unweighted-list"
            else create_weighted_loading_df(genecart_id)
        )
    except Exception as e:
        print(str(e), file=fh)
        traceback.print_exc()
        status["status"] = "failed"
        status["error"] = str(e)
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    # Assumes first column is unique identifiers. Standardize on a common index name
    loading_df = loading_df.rename(columns={loading_df.columns[0]: "dataRowNames"})

    # Drop the gene symbol column
    loading_df = loading_df.drop(loading_df.columns[1], axis=1)

    loading_df = loading_df.set_index("dataRowNames")

    # If cross-species, remap the genecart genes to the orthologous genes for the dataset's organism
    try:
        # Get the organism ID associated with the dataset.
        ds = geardb.get_dataset_by_id(dataset_id)
        if not ds:
            raise Exception(
                "Dataset not found in database. Please contact a gEAR admin."
            )

        # Both the genecart and dataset need to have organism IDs, which should have happened at upload time
        if not genecart.organism_id:
            raise Exception(
                "Gene cart {} does not have an organism ID. Please contact a gEAR admin.".format(
                    genecart_id
                )
            )

        if not ds.organism_id:
            raise Exception(
                "Dataset {} does not have an organism ID. Please contact a gEAR admin.".format(
                    dataset_id
                )
            )

        if not genecart.organism_id == ds.organism_id:
            ortholog_file = get_ortholog_file(
                str(genecart.organism_id), str(ds.organism_id), ANNOTATION_TYPE
            )
            if ortholog_file is None:
                raise Exception(
                    "Could not find an orthologous mapping file between the gene list organism and the dataset organism."
                )
            loading_df = map_dataframe_genes(loading_df, ortholog_file)
    except Exception as e:
        print(str(e), file=fh)
        traceback.print_exc()
        status["status"] = "failed"
        status["success"] = -1
        status["error"] = str(e)
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    # Drop duplicate unique identifiers. This may happen if two unweighted gene cart genes point to the same Ensembl ID in the db
    loading_df = loading_df[~loading_df.index.duplicated(keep="first")]

    is_spatial = ds.dtype == "spatial"

    # NOTE Currently no analyses are supported yet.
    try:
        ana = get_analysis(None, dataset_id, session_id, is_spatial)
    except Exception:
        traceback.print_exc()
        status["status"] = "failed"
        status["error"] = "Could not retrieve analysis."
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    try:
            args = {}
            if is_spatial:
                args['include_images'] = False
            else:
                args['backed'] = True
            adata = ana.get_adata(**args)
    except Exception:
        traceback.print_exc()
        status["status"] = "failed"
        status["error"] = "Could not retrieve analysis data."
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    # If dataset genes have duplicated index names, we need to rename them to avoid errors
    # in collecting rownames in projectR (which gives invalid output)
    # This means these duplicated genes will not be in the intersection of the dataset and pattern genes
    if isinstance(ana, SpatialAnalysis):
        dedup_copy = str(ana.dataset_path().replace(".zarr", ".dups_removed.h5ad"))
    else:
        dedup_copy = str(ana.dataset_path().replace(".h5ad", ".dups_removed.h5ad"))
    dedup_copy = Path(dedup_copy)

    if (adata.var.index.duplicated(keep="first")).any():
        if dedup_copy.exists():
            dedup_copy.unlink()
        adata = adata[:, ~adata.var.index.duplicated(keep="first")].copy(
            filename=dedup_copy
        )

    num_target_genes = adata.shape[1]
    num_loading_genes = loading_df.shape[0]

    # Perform overlap to see if there are overlaps between genes from both dataframes
    index_intersection = adata.var.index.intersection(loading_df.index)
    intersection_size = index_intersection.size

    if intersection_size == 0:
        message = "No common genes between the target dataset ({} genes) and the pattern file ({} genes).".format(
            num_target_genes, num_loading_genes
        )
        status["status"] = "failed"
        status["error"] = message
        status["result"] =  {
            "success": -1,
            "num_common_genes": intersection_size,
            "num_genecart_genes": num_loading_genes,
            "num_dataset_genes": num_target_genes,
        }
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    message = "Found {} common genes between the target dataset ({} genes) and the pattern file ({} genes).".format(
        intersection_size, num_target_genes, num_loading_genes
    )

    # Reduce the size of both dataframes before POSTing
    adata = adata[:, index_intersection]  # type: ignore
    loading_df = loading_df.loc[index_intersection]
    loading_df = loading_df.fillna(0)  # Fill NaN values with 0

    # Ensure target dataset has genes as rows
    # ! The to_df() seems to be a memory_chokepoint, doubling the memory used by the adata object
    target_df = adata.to_df().transpose()

    # If zscore is enabled, scale the target_df according to zscore by row (genes)
    # For NaN values, they are ignored in the calculation
    if zscore:
        target_df = target_df.apply(
            lambda row: pd.Series(stats.zscore(row, ddof=1, nan_policy="omit"), index=row.index),
            axis=1
        )

    target_df = target_df.fillna(0)  # Fill NaN values with 0

    # Close dataset adata so that we do not have a stale opened object
    adata.file.close()

    dataset_projection_csv = build_projection_csv_path(
        dataset_id, projection_id, "dataset"
    )

    # This is about the time I could consider this to be running, as it is the start of the "long-running" part of the task
    status["status"] = "running"
    write_projection_status(JOB_STATUS_FILE, status)

    # Create lock file if it does not exist
    lockfile = str(dataset_projection_csv) + ".lock"
    if Path(lockfile).exists():
        print(
            "INFO: Found lockfile for another current projectR run of {}.  Going to wait for that run to finish and steal its output.".format(
                projection_id
            ),
            file=fh,
        )
        try:
            # Test to see if the exclusive lock has expired
            lock_fh = create_lock_file(lockfile)
            print("INFO: Lock for {} seems to be stale. Removing it.".format(projection_id), file=fh)
            remove_lock_file(lock_fh, lockfile)
        except Exception:
            print("INFO: Lock for {} seems to be valid.".format(projection_id), file=fh)
            # If lock belongs to a valid run, wait for lock to be removed,
            # then return info that is normally returned after projectR is run
            while True:
                sleep(1)
                if not Path(lockfile).exists():
                    status["result"] = {
                        "success": 2,
                        "message": message,
                        "projection_id": projection_id,
                        "num_common_genes": intersection_size,
                        "num_genecart_genes": num_loading_genes,
                        "num_dataset_genes": num_target_genes,
                    }
                    return status

    try:
        lock_fh = create_lock_file(lockfile)
        # NOTE: Will not trigger if process crashes or is forcibly killed off.
    except IOError:
        # This should ideally never be encountered as the previous code should handle existing locked files
        message = "Could not create lock file for this projectR run."
        status["status"] = "failed"
        status["error"] = message
        status["result"] = {
            "success": -1,
            "num_common_genes": intersection_size,
            "num_genecart_genes": num_loading_genes,
            "num_dataset_genes": num_target_genes,
        }
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    # Chunk size needs to adjusted by how many genes are present, so that the payload always stays under the body size limit
    chunk_size = calculate_chunk_size(len(target_df.index), len(target_df.columns))
    if len(target_df.columns) < chunk_size:
        chunk_size = len(target_df.columns)

    print(
        "TARGET: {}\nGENECART: {}\nTARGET DF (genes,samples): {}\nSAMPLES PER CHUNK: {}".format(
            dataset_id, genecart_id, target_df.shape, chunk_size
        ),
        file=fh,
    )

    # report number of chunks to make
    index_slices = sliced(range(len(target_df.columns)), chunk_size)
    print("NUMBER OF CHUNKS: {}".format(len(list(index_slices))), file=fh)

    # shuffle target dataframe rows.  Needed to balance out the chunks in the NMF algorithms. Seeded for reproducibility.
    target_df = target_df.sample(frac=1, random_state=42)

    if algorithm == "fixednmf":
        # Normalize the target_df by the minimum sum of each expression column
        columns_sums = target_df.sum(axis=0)
        normalized_columns_sums = columns_sums / columns_sums.min()
        target_df = target_df.div(normalized_columns_sums, axis=1)

    projection_pval_df = pd.DataFrame()

    if this.servercfg["projectR_service"]["cloud_run_enabled"].startswith("1"):
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        try:
            loop.run_until_complete(
                fetch_all_queue(
                    target_df,
                    loading_df,
                    algorithm,
                    full_output,
                    projection_id,
                    chunk_size,
                    fh,
                    concurrency=CONCURRENT_REQUEST_LIMIT,
                )
            )
            print("INFO: All fetch tasks have completed", file=fh)
        except asyncio.TimeoutError:
            remove_lock_file(lock_fh, lockfile)
            status["status"] = "failed"
            status["error"] = "Timeout while waiting for projectR to complete."
            status["result"] = {
                "success": -1,
                "num_common_genes": intersection_size,
                "num_genecart_genes": num_loading_genes,
                "num_dataset_genes": num_target_genes,
            }
            return status
        except Exception as e:
            print(str(e), file=fh)
            # Raises as soon as one "gather" task has an exception
            remove_lock_file(lock_fh, lockfile)
            status["status"] = "failed"
            status["error"] = "Something went wrong with the projection-creating step."
            status["result"] = {
                "success": -1,
                "num_common_genes": intersection_size,
                "num_genecart_genes": num_loading_genes,
                "num_dataset_genes": num_target_genes,
            }
            return status
        finally:
            write_projection_status(JOB_STATUS_FILE, status)
            # Wait 250 ms for the underlying SSL connections to close
            loop.run_until_complete(asyncio.sleep(0.250))
            loop.stop()  # prevent "Task was destroyed but it is pending!" messages
            loop.close()

        print("INFO: Concatenating results to dataframe", file=fh)

        # single result = {"projection": "json", "pval": "json"}
        # concatenate all the results into a single DataFrame for each key

        projection_patterns_df, projection_pval_df = concat_fetch_results_to_dataframe(
            projection_id
        )

        gc.collect()  # trying to clear memory

        if len(projection_patterns_df.index) != len(adata.obs.index):
            message = "Not all chunked sample rows were returned by projectR. Saved partial results to disk. Refresh to try again."
            print(message, file=fh)
            remove_lock_file(lock_fh, lockfile)
            status["status"] = "failed"
            status["error"] = message
            status["result"] = {
                "success": -1,
                "num_common_genes": intersection_size,
                "num_genecart_genes": num_loading_genes,
                "num_dataset_genes": num_target_genes,
            }
            write_projection_status(JOB_STATUS_FILE, status)
            return status

        # There is a good chance the samples are now out of order, which will break
        # the copying of the dataset observation metadata when this output is converted
        # to an AnnData object. So reorder back to dataset sample order.
        projection_patterns_df = projection_patterns_df.reindex(
            adata.obs.index.tolist()
        )
        if not projection_pval_df.empty:
            projection_pval_df = projection_pval_df.reindex(adata.obs.index.tolist())

        # Delete all the chunk output files for this projection_id
        for filepath in CHUNK_OUTPUTS_DIR.glob(f"{projection_id}_chunk*.json"):
            try:
                filepath.unlink()
            except Exception as e:
                print(
                    "WARNING: Could not delete chunk output file {}: {}".format(
                        filepath, str(e)
                    ),
                    file=fh,
                )

    else:
        # If not using the cloud run service, do this on the server
        abs_path_gear = Path(__file__).resolve().parents[3]
        service_path = abs_path_gear.joinpath("services")
        # abs_path_lib is a Path object so we need to convert to string
        sys.path.insert(0, str(service_path))

        # https://github.com/IGS/gEAR/issues/442#issuecomment-1317239909
        # Basically this is a stopgap until projectR has an option to remove
        # centering around zero for PCA loadings.  Chunking the data breaks
        # the output due to the centering around zero step.
        try:
            if algorithm == "pca":
                from projectr.main import do_pca_projection

                projection_patterns_df = do_pca_projection(target_df, loading_df)
            elif algorithm == "binary":
                from projectr.main import do_binary_projection

                projection_patterns_df = do_binary_projection(target_df, loading_df)
            elif algorithm == "2silca":
                pass
            elif algorithm in ["nmf", "fixednmf"]:
                from projectr.rfuncs import run_projectR_cmd

                # R code: projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix)
                # "projection" is a DataFrame (index 0)
                # "pval" is a matrix (index 1)

                projection_patterns = run_projectR_cmd(
                    target_df, loading_df, algorithm, full_output
                )

                projection_patterns_df = projection_patterns[0].transpose()

                if full_output and algorithm == "nmf":
                    projection_pval_df = projection_patterns[1].transpose()

                projection_patterns = run_projectR_cmd(
                    target_df, loading_df, algorithm, full_output
                )

            else:
                raise ValueError("Algorithm {} is not supported".format(algorithm))
        except Exception as e:
            # clear lock file
            remove_lock_file(lock_fh, lockfile)
            print(str(e), file=fh)
            status["status"] = "failed"
            status["error"] = "Something went wrong with the projection-creating step."
            status["result"] = {
                "success": -1,
                "num_common_genes": intersection_size,
                "num_genecart_genes": num_loading_genes,
                "num_dataset_genes": num_target_genes,
            }
            return status
        finally:
            write_projection_status(JOB_STATUS_FILE, status)

    # Close adata so that we do not have a stale opened object
    if adata.isbacked:
        adata.file.close()

    if dedup_copy.exists():
        dedup_copy.unlink()

    # Have had cases where the column names are x1, x2, x3, etc. so load in the original pattern names
    projection_patterns_df = projection_patterns_df.set_axis(
        loading_df.columns, axis="columns"
    )
    if not projection_pval_df.empty:
        projection_pval_df = projection_pval_df.set_axis(
            loading_df.columns, axis="columns"
        )

    # Check that all DataFrame values are not null.  If not, we cannot proceed.
    # Ultimately, we cannot plot this and do not want to write to file.
    # This could be due to an R code error or something in post-processing
    if projection_patterns_df.isna().to_numpy().any():
        message = "There are NaN values in the projection patterns.  Cannot proceed."
        remove_lock_file(lock_fh, lockfile)
        status["status"] = "failed"
        status["error"] = message
        status["result"] = {
            "success": -1,
            "num_common_genes": intersection_size,
            "num_genecart_genes": num_loading_genes,
            "num_dataset_genes": num_target_genes,
        }
        write_projection_status(JOB_STATUS_FILE, status)
        return status

    # if full_output = True, then write the pval matrix to file using the same name as the projection file
    if full_output and algorithm == "nmf":
        if projection_pval_df.empty:
            status["status"] = "failed"
            status["error"] = "No pval matrix was generated by projectR."
            status["result"] = {
                "success": -1,
                "num_common_genes": intersection_size,
                "num_genecart_genes": num_loading_genes,
                "num_dataset_genes": num_target_genes,
            }
            write_projection_status(JOB_STATUS_FILE, status)
            return status
        projection_pval_csv = build_projection_csv_path(
            dataset_id, projection_id, "pval"
        )
        projection_pval_df.to_csv(projection_pval_csv)

    print(
        "INFO: Writing projection patterns to {}".format(dataset_projection_csv),
        file=fh,
    )
    projection_patterns_df.to_csv(dataset_projection_csv)

    del projection_patterns_df
    gc.collect()  # trying to clear memory

    dataset_projection_json_file = build_projection_json_path(dataset_id, "dataset")

    # Add new configuration to the list for this dictionary key
    with open(dataset_projection_json_file) as projection_fh:
        try:
            dataset_projections_dict = json.load(projection_fh)
        except json.JSONDecodeError:
            dataset_projections_dict = {}

    dataset_projections_dict.setdefault(genecart_id, []).append(
        {
            "uuid": projection_id,
            "algorithm": algorithm,
            "zscore": zscore,
            "num_common_genes": intersection_size,
            "num_genecart_genes": num_loading_genes,
            "num_dataset_genes": num_target_genes,
        }
    )
    write_to_json(dataset_projections_dict, dataset_projection_json_file)

    # Do the same for the genecart version
    genecart_projection_csv = build_projection_csv_path(
        genecart_id, projection_id, "genecart"
    )
    genecart_projection_json_file = build_projection_json_path(genecart_id, "genecart")

    # Symlink dataset_projection_csv to genecart_projection_csv (this syntax feels like it's in reverse but it is correct)
    # NOTE: In the Docker instance, symlink reflects path to the mounted volume, not the local path
    try:
        genecart_projection_csv.symlink_to(dataset_projection_csv)
    except FileExistsError:
        print("Symlink already exists for {}".format(dataset_projection_csv), file=fh)

    with open(genecart_projection_json_file) as projection_fh:
        try:
            genecart_projections_dict = json.load(projection_fh)
        except json.JSONDecodeError:
            genecart_projections_dict = {}
    genecart_projections_dict.setdefault(dataset_id, []).append(
        {
            "uuid": projection_id,
            "algorithm": algorithm,
            "zscore": zscore,
            "num_common_genes": intersection_size,
            "num_genecart_genes": num_loading_genes,
            "num_dataset_genes": num_target_genes,
        }
    )
    write_to_json(genecart_projections_dict, genecart_projection_json_file)

    # Remove file lock
    print("INFO: Removing lock file for {}".format(projection_id), file=fh)
    remove_lock_file(lock_fh, lockfile)

    status["status"] = "complete"
    status["result"] = {
        "success": success,
        "message": message,
        "projection_id": projection_id,
        "num_common_genes": intersection_size,
        "num_genecart_genes": num_loading_genes,
        "num_dataset_genes": num_target_genes,
    }
    write_projection_status(JOB_STATUS_FILE, status)
    return status


class ProjectROutputFile(Resource):
    """
    Get or create path to output file. Will mkdir if the directory does not currently exist.
    """

    def post(self, dataset_id: str) -> dict:
        args = parser.parse_args()
        genecart_id = args["genecart_id"]
        algorithm = args["algorithm"]
        zscore = args["zscore"]

        if not zscore:
            zscore = False

        # Create the directory if it doesn't exist
        # NOTE: The "mkdir" and "touch" commands are not atomic. Do not fail if directory or file was created by another process.
        dataset_projection_json_file = build_projection_json_path(dataset_id, "dataset")
        if not dataset_projection_json_file.parent.is_dir():
            dataset_projection_json_file.parent.mkdir(parents=True, exist_ok=True)

        genecart_projection_json_file = build_projection_json_path(
            genecart_id, "genecart"
        )
        if not genecart_projection_json_file.parent.is_dir():
            genecart_projection_json_file.parent.mkdir(parents=True, exist_ok=True)

        # Create the genecart projection json file if it doesn't exist
        # We will use the dataset json file as a check if a projection already exists though
        if not genecart_projection_json_file.is_file():
            genecart_projection_json_file.touch(exist_ok=True)
            write_to_json({}, genecart_projection_json_file)  # create empty json file

        # If the dataset json file exists, we can read it and get the output file path
        if not dataset_projection_json_file.is_file():
            dataset_projection_json_file.touch(exist_ok=True)
            write_to_json({}, dataset_projection_json_file)  # create empty json file
            return {"projection_id": None}

        with open(dataset_projection_json_file) as projection_fh:
            projections_dict = json.load(projection_fh)

        # If the pattern was not projected onto this dataset, initialize a list of configs
        # "cart.{projection_id}" added for backwards compatability
        if not (
            genecart_id in projections_dict
            or "cart.{}".format(genecart_id) in projections_dict
        ):
            # NOTE: tried to write JSON with empty list but it seems that empty keys are skipped over.
            return {"projection_id": None}

        # If legacy version exists, copy to current format.
        if (genecart_id not in projections_dict) and "cart.{}".format(
            genecart_id
        ) in projections_dict:
            print(
                "Copying legacy cart.{} to {} in the projection json file.".format(
                    genecart_id, genecart_id
                ),
                file=sys.stderr,
            )
            projections_dict.pop("cart.{}".format(genecart_id), None)
            # move legacy genecart projection stuff to new version
            print(
                "Moving legacy cart.{} contents to {} in the projection genecart directory.".format(
                    genecart_id, genecart_id
                ),
                file=sys.stderr,
            )
            old_genecart_projection_json_file = build_projection_json_path(
                "cart.{}".format(genecart_id), "genecart"
            )
            old_genecart_projection_json_file.parent.rename(
                dataset_projection_json_file.parent
            )

        for config in projections_dict[genecart_id]:
            # Handle legacy JSON outputs (i.e. old or non-existing names)
            if "is_pca" in config:
                config["algorithm"] = "pca" if config["is_pca"] == 1 else "nmf"
            if "zscore" not in config:
                config["zscore"] = False

            # Check if the config matches the current request (full output does not matter here)
            if algorithm == config["algorithm"] and zscore == config["zscore"]:
                projection_id = config["uuid"]
                dataset_projection_csv = build_projection_csv_path(
                    dataset_id, projection_id, "dataset"
                )
                return {
                    "projection_id": projection_id
                    if Path(dataset_projection_csv).is_file()
                    else None
                }

        # If we get here, we didn't find a file for this config
        return {"projection_id": None}


class ProjectR(Resource):
    """
    Formats inputs to prep for projectR API call running on Google Cloud Run
    """

    def post(self, dataset_id: str) -> dict:
        session_id = request.cookies.get("gear_session_id", "")
        args = run_projectr_parser.parse_args()

        genecart_id = args["genecart_id"]
        algorithm = args["algorithm"]
        projection_id = args["projection_id"]
        scope = args["scope"]
        zscore = args["zscore"]
        full_output = args["full_output"]

        # Use default based on gear.ini settings
        # This should be default to True coming from the UI but this would be a nice fallback.
        if full_output is None:
            if "full_output" in this.servercfg["projectR_service"]:
                full_output = this.servercfg["projectR_service"]["full_output"].startswith("1")

        # Currently only NMF runs through the actual projectR code and can give full output
        if algorithm not in ["nmf", "fixednmf"]:
            full_output = False

        uuid_args = (dataset_id, genecart_id, algorithm, zscore)

        projection_id = projection_id or str(create_new_uuid(*uuid_args))
        dataset_projection_csv = build_projection_csv_path(
            dataset_id, projection_id, "dataset"
        )
        dataset_projection_json_file = build_projection_json_path(dataset_id, "dataset")

        run_projectr = True

        status = init_job_status(projection_id)

        # Housekeeping... create some dir paths if they do not exist
        JOB_STATUS_DIR.mkdir(parents=True, exist_ok=True)
        CHUNK_OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)

        # If projectR has already been run, we can just load the csv file.  Otherwise, let it rip!
        if Path(dataset_projection_csv).is_file():
            print(
                "INFO: Found exisitng dataset_projection_csv file {}, loading it.".format(
                    dataset_projection_csv
                )
            )

            run_projectr = False
            if full_output:
                # does pval matrix exist?
                pval_csv = build_projection_csv_path(dataset_id, projection_id, "pval")
                if not Path(pval_csv).is_file():
                    print(
                        "INFO: Full output requested but pval matrix file {} does not exist.".format(
                            pval_csv
                        )
                    )
                    run_projectr = True

            if not run_projectr:
                # Projection already exists, so we can just return info we want to return in a message
                with open(dataset_projection_json_file) as projection_fh:
                    projections_dict = json.load(projection_fh)
                common_genes = None
                genecart_genes = None
                dataset_genes = None

                # Get projection info from the json file if it already exists
                if genecart_id in projections_dict:
                    for config in projections_dict[genecart_id]:
                        # Handle legacy JSON outputs (i.e. old or non-existing names)
                        if "is_pca" in config:
                            config["algorithm"] = (
                                "pca" if config["is_pca"] == 1 else "nmf"
                            )
                        if "zscore" not in config:
                            config["zscore"] = False

                        # Check if the config matches the current request
                        if (
                            algorithm == config["algorithm"]
                            and zscore == config["zscore"]
                        ):
                            common_genes = config.get("num_common_genes", None)
                            genecart_genes = config.get("num_genecart_genes", -1)
                            dataset_genes = config.get("num_dataset_genes", -1)
                            break

                message = ""
                if common_genes:
                    message = "Found {} common genes between the target dataset ({} genes) and the pattern file ({} genes).".format(
                        common_genes, dataset_genes, genecart_genes
                    )

                result = {
                    "success": 1,
                    "message": message,
                    "projection_id": projection_id,
                    "num_common_genes": common_genes,
                    "num_genecart_genes": genecart_genes,
                    "num_dataset_genes": dataset_genes,
                }
                status["status"] = "complete"
                status["result"] = result
                return status

        JOB_STATUS_FILE = JOB_STATUS_DIR.joinpath(f"job_{projection_id}.json")
        # if this file does not exist, write the status
        # If it does exist, use this projections file's run as our own.
        if Path(JOB_STATUS_FILE).is_file():
            with open(JOB_STATUS_FILE, "r") as fh:
                status = json.load(fh)
                if status["status"] in ["pending", "running", "complete"]:
                    print(f"[x] Job {projection_id} is already {status['status']}", file=sys.stderr)
                    return status
                elif status["status"] == "failed":
                    # delete status file so we can start a rerun
                    print(f"[x] Job {projection_id} has failed. Attempting a rerun", file=sys.stderr)
                    Path(JOB_STATUS_FILE).unlink(missing_ok=True)
                    # Ensure "error" status is not written to file for new polling session
                    status = init_job_status(projection_id)

        # Write pending state
        write_projection_status(JOB_STATUS_FILE, status)

        # Create a messaging queue if necessary. Make it persistent across the lifetime of the Flask server.
        # Channels will be spawned during each task.
        if this.servercfg["projectR_service"]["queue_enabled"].startswith("1"):

            import gearqueue
            host = this.servercfg["projectR_service"]["queue_host"]

            try:
                # Connect as a blocking RabbitMQ publisher
                connection = gearqueue.Connection(
                    host=host, publisher_or_consumer="publisher"
                )
            except Exception as e:
                status["status"] = "failed"
                status["error"] = str(e)
                traceback.print_exc(file=sys.stderr)
                return status

            # Connect as a blocking RabbitMQ publisher
            with connection:
                connection.open_channel()

                # Create the publisher
                payload = dict()
                payload["dataset_id"] = dataset_id
                payload["session_id"] = session_id
                payload["projection_id"] = (
                    projection_id  # the "arg" was assigned and potentially modified
                )
                payload["genecart_id"] = genecart_id
                payload["scope"] = scope
                payload["algorithm"] = algorithm
                payload["zscore"] = zscore
                payload["full_output"] = full_output

                try:
                    connection.publish(
                        queue_name="projectr",
                        message=payload,  # method dumps JSON
                    )
                    print(
                        "[x] Requesting for dataset {} and genecart {}".format(
                            dataset_id, genecart_id
                        ),
                        file=sys.stderr,
                    )

                    status["result"] = {"projection_id": projection_id}
                    return status
                except Exception as e:
                    status["status"] = "failed"
                    status["error"] = str(e)
                    traceback.print_exc(file=sys.stderr)
                    write_projection_status(JOB_STATUS_FILE, status)
                    return status

        else:
            status = projectr_callback(
                dataset_id,
                genecart_id,
                projection_id,
                session_id,
                scope,
                algorithm,
                zscore,
                full_output,
                sys.stderr,
            )
            # Delete job status file
            Path(JOB_STATUS_FILE).unlink(missing_ok=True)
            return status


class ProjectRStatus(Resource):
    """
    Get the status of a ProjectR job.
    """
    def get(self, projection_id):
        safe_projection_id = secure_filename(str(projection_id))
        JOB_STATUS_FILE = JOB_STATUS_DIR.joinpath(f"job_{safe_projection_id}.json")
        # Validate the final path is within the job status dir
        resolved_status_file = JOB_STATUS_FILE.resolve()
        if not resolved_status_file.is_relative_to(JOB_STATUS_DIR.resolve()):
            # Reject attempts to escape the directory
            from flask import abort
            abort(403, description="Invalid job id/path")
        with open(resolved_status_file, "r") as fh:
            status = json.load(fh)

        # possible states - running, complete, failed
        if status["status"] == "complete":
            # Delete job status file
            Path(JOB_STATUS_FILE).unlink(missing_ok=True)

        return status

