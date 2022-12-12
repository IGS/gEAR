#!/opt/bin/python3

"""
This script is to process a single projectR-to-tsne plot run in order to do some memory profiling.

Shaun Adkins
"""

import json, sys, fcntl
import pandas as pd

from os import getpid
from time import sleep
from more_itertools import sliced

import asyncio


import base64
import io
import os
import re
from math import ceil
from pathlib import Path

lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)
import geardb

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from matplotlib import cm

from memory_profiler import profile

sc.settings.set_figure_params(dpi=100)
sc.settings.verbosity = 0
sc.settings.figdir = '/tmp/'

# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

dataset_id = "1b12dde9-1762-7564-8fbd-1b07b750505f"
genecart_id = "0a0216f1"
is_pca = True
session_id = 'ee95e48d-c512-4083-bf05-ca9f65e2c12a'
analysis_owner_id = "662"
gene_symbol = "PC1"
plot_type = "tsne_static"
x_axis = "tSNE_1"
y_axis = "tSNE_2"
colorize_by = "cell_type"

#dataset_id = "8dbfdcd7-a826-4658-e73a-bf85fabe7d6b"
#genecart_id = "79189c08"
#is_pca = False
#session_id = '5c3898f5-85ea-4fdc-b97a-496a93d489d6'
#analysis_owner_id = "483"
#gene_symbol = "p1"
#plot_type = "umap_static"
#x_axis = "UMAP_0"
#y_axis = "UMAP_1"
#colorize_by = "cell.type"


scope = "weighted-list"
user = geardb.get_user_from_session_id(session_id)
analysis = None
colorblind_mode = False
colors = {}
horizontal_legend = False
max_columns = None
order = {}
plot_by_group = None
projection_id = None
skip_gene_plot = None

ONE_LEVEL_UP = 1
abs_path_git = Path(__file__).resolve().parents[ONE_LEVEL_UP] # git root
abs_path_www = abs_path_git.joinpath("www")
CARTS_BASE_DIR = abs_path_www.joinpath("carts")
PROJECTIONS_BASE_DIR = abs_path_www.joinpath("projections")
ORTHOLOG_BASE_DIR = abs_path_www.joinpath("feature_mapping")
PROJECTIONS_JSON_BASENAME = "projections.json"

ANNOTATION_TYPE = "ensembl" # NOTE: This will change in the future to be varied.

# limit of asynchronous tasks that can happen at a time
# I am setting this slightly under the "MaxKeepAliveRequests" in apache.conf
SEMAPHORE_LIMIT = 50

PLOT_TYPE_TO_BASIS = {
    "tsne_static": "tsne",
    "tsne": "tsne",             # legacy
    "umap_static": "umap",
    "pca_static": "pca"
}
COLOR_HEX_PTRN = r"^#(?:[0-9a-fA-F]{3}){1,2}$"

NUM_LEGENDS_PER_COL = 12    # Max number of legend items per column allowed in vertical legend
NUM_HORIZONTAL_COLS = 8 # Number of columns in horizontal legend

class PlotError(Exception):
    """Error based on plotting issues."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)


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
    res_dfs = [pd.read_json(res_json, orient="split", dtype="float32") for res_json in res_jsons]
    projection_patterns_df = pd.concat(res_dfs)
    return projection_patterns_df

def create_new_uuid(dataset_id, genecart_id, is_pca):
    import hashlib, uuid
    uuid_str = "{}-{}-{}".format(dataset_id, genecart_id, is_pca)
    md5 = hashlib.md5()
    md5.update(uuid_str.encode("utf-8"))
    return uuid.UUID(md5.hexdigest())

def create_unweighted_loading_df(gc):
    # Now convert into a GeneCollection to get the Ensembl IDs (which will be the unique identifiers)
    gene_collection = geardb.GeneCollection()
    gene_collection.get_by_gene_symbol(gene_symbol=" ".join(gc.genes), exact=True)

    loading_data = []
    for gene in gene_collection.genes:
        loading_data.append({"dataRowNames": gene.ensembl_id, "gene_sym":gene.gene_symbol, "unweighted":1})
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

async def fetch_all(loop, target_df, loading_df, is_pca, genecart_id, dataset_id, chunk_size):
    """Create coroutine tasks out of all chunked projection cloud run service POST requests."""
    # Code influenced by https://stackoverflow.com/a/68288374 and https://superfastpython.com/asyncio-as_completed/
    import aiohttp
    sem = asyncio.Semaphore(SEMAPHORE_LIMIT) # limit simultaneous tasks so the gEAR server CPU isn't overloaded
    async with aiohttp.ClientSession(loop=loop) as client, sem:
        # Create coroutines to be executed.
        # Wrap in "asyncio.create_task" to run concurrently.
        coros = [asyncio.create_task(fetch_one(client, {
                "target": chunk_df.to_json(orient="split")
                , "loadings": loading_df.to_json(orient="split")
                , "is_pca": is_pca
                , "genecart_id":genecart_id # This helps in identifying which combinations are going through
                , "dataset_id":dataset_id
                })) for chunk_df in chunk_dataframe(target_df, chunk_size)]
        # This loop processes results as they come in.
        # asyncio.as_completed creates a generator from the coroutines/tasks
        return [await coro for coro in asyncio.as_completed(coros)]

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
        return await response.json()

def calculate_figure_height(num_plots):
    """Determine height of tsne plot based on number of group elements."""
    return (num_plots * 4) + (num_plots -1)

def calculate_figure_width(num_plots):
    """Determine width of tsne plot based on number of group elements."""
    return (num_plots * 6) + (num_plots -1)

def calculate_num_legend_cols(group_len):
    """Determine number of columns legend should have in tSNE plot."""
    return ceil(group_len / NUM_LEGENDS_PER_COL)

def create_colorscale_with_zero_gray(colorscale):
    """Take a predefined colorscale, and change the 0-value color to gray, and return."""
    from matplotlib import cm
    from matplotlib.colors import ListedColormap

    # Create custom colorscale with gray at the 0.0 level
    # Src: https://matplotlib.org/tutorials/colors/colormap-manipulation.html
    ylorrd = cm.get_cmap(colorscale, 256)
    newcolors = ylorrd(np.linspace(0, 1, 256))  # split colormap into 256 parts over 0:1 range
    gray = np.array([192/256, 192/256, 192/256, 1])
    newcolors[0, :] = gray
    return mcolors.ListedColormap(newcolors)

def create_projection_adata(dataset_adata, dataset_id, projection_id):
    # Create AnnData object out of readable CSV file
    # ? Does it make sense to put this in the geardb/Analysis class?
    projection_dir = Path(PROJECTIONS_BASE_DIR).joinpath("by_dataset", dataset_id)
    projection_adata_path = projection_dir.joinpath("{}.h5ad".format(projection_id))
    if projection_adata_path.is_file():
        os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'   # see https://github.com/h5py/h5py/issues/1722
        return sc.read_h5ad(projection_adata_path, backed="r")

    projection_csv_path = projection_dir.joinpath("{}.csv".format(projection_id))
    try:
        projection_adata = sc.read_csv(projection_csv_path)
    except:
        raise PlotError("Could not create projection AnnData object from CSV.")
    projection_adata.obs = dataset_adata.obs
    # Close dataset adata so that we do not have a stale opened object
    if dataset_adata.isbacked:
        dataset_adata.file.close()
    projection_adata.var["gene_symbol"] = projection_adata.var_names
    # Associate with a filename to ensure AnnData is read in "backed" mode
    projection_adata.filename = projection_adata_path
    return projection_adata

def get_colorblind_scale(n_colors):
    """Get a colorblind friendly colorscale (Viridis). Return n colors spaced equidistantly."""
    cividis = cm.get_cmap("viridis", n_colors)
    colors = cividis.colors
    # convert to hex since I ran into some issues using rpg colors
    return [mcolors.rgb2hex(color) for color in colors]

def sort_legend(figure, sort_order, horizontal_legend=False):
    """Sort legend of plot."""
    handles, labels = figure.get_legend_handles_labels()
    new_handles = [handles[idx] for idx, name in enumerate(sort_order)]
    new_labels = [labels[idx] for idx, name in enumerate(sort_order)]

    # If horizontal legend, we need to sort in a way to have labels read from left to right
    if horizontal_legend:
        leftover_cards = len(new_handles) % NUM_HORIZONTAL_COLS
        num_chunks = int(len(new_handles) / NUM_HORIZONTAL_COLS)

        # If number of groups is less than num_cols, they can just be put on a single line
        if num_chunks == 0:
            return (new_handles, new_labels)

        # Split into relatively equal chumks.
        handles_sublists = np.array_split(np.array(new_handles), num_chunks)
        labels_sublists = np.array_split(np.array(new_labels), num_chunks)

        # Zipping gets weird if there's a remainder so remove those leftover items and add back later
        if leftover_cards:
            handles_leftover = new_handles[-leftover_cards:]
            labels_leftover = new_labels[-leftover_cards:]

            handles_sublists = np.array_split(np.array(new_handles[:-leftover_cards]), num_chunks)
            labels_sublists = np.array_split(np.array(new_labels[:-leftover_cards]), num_chunks)

        # Zip numpy arrays, then flatten into a 1D list.  Add leftover cards as well to end
        new_handles = np.column_stack(handles_sublists).flatten().tolist()
        new_labels = np.column_stack(labels_sublists).flatten().tolist()

        # Insert leftover cards back into the list. Start from back to front so indexes are not manipulated.
        for i in reversed(range(leftover_cards)):
            new_handles.insert(num_chunks * (i+1), handles_leftover[i])
            new_labels.insert(num_chunks * (i+1), labels_leftover[i])

    return (new_handles, new_labels)

def run_projection_prestep():
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

@profile
def run_projection():
    success = 1
    message = ""

    global projection_id
    projection_id = projection_id or create_new_uuid(dataset_id, genecart_id, is_pca)
    dataset_projection_csv = build_projection_csv_path(dataset_id, projection_id, "dataset")
    dataset_projection_json_file = build_projection_json_path(dataset_id, "dataset")

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
            message = "Found {} common genes between the target dataset ({} genes) and the pattern file ({} genes).".format(common_genes, dataset_genes, genecart_genes)

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
        * Output - Rows are patterns, Cols are data observations (before transposing)

    Only step 3 needs to be in R, and we use rpy2 to call that.
    """

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

    # Row: Genes
    # Col: Pattern weights
    try:
        loading_df = create_unweighted_loading_df(gc) if scope == "unweighted-list" else create_weighted_loading_df(genecart_id)
    except Exception as e:
        print(str(e), file=sys.stderr)
        return {
            'success': -1
            , 'message': str(e)
        }

    # Assumes first column is unique identifiers. Standardize on a common index name
    loading_df.rename(columns={ loading_df.columns[0]:"dataRowNames" }, inplace=True)

    # Drop the gene symbol column
    loading_df.drop(loading_df.columns[1], axis=1, inplace=True)

    loading_df.set_index('dataRowNames', inplace=True)

    # If cross-species, remap the genecart genes to the orthologous genes for the dataset's organism
    try:
        # Get the organism ID associated with the dataset.
        ds = geardb.get_dataset_by_id(dataset_id)
        if not gc.organism_id == ds.organism_id:
            orthomap_file_base = "orthomap.{0}.{2}__{1}.{2}.hdf5".format(gc.organism_id, ds.organism_id, ANNOTATION_TYPE)
            orthomap_file = ORTHOLOG_BASE_DIR.joinpath(orthomap_file_base)
            try:
                loading_df = remap_df_genes(loading_df, orthomap_file)
            except:
                message = "Could not remap pattern genes to ortholog equivalent in the dataset"
                return {
                    "success": -1
                    , "message": message
                }
    except:
        message = "Dataset was not found in the database. Please contact a gEAR admin."
        return {"success": -1, "message": message}

    # Drop duplicate unique identifiers. This may happen if two unweighted gene cart genes point to the same Ensembl ID in the db
    loading_df = loading_df[~loading_df.index.duplicated(keep='first')]

    # NOTE Currently no analyses are supported yet.
    try:
        ana = get_analysis(None, dataset_id, session_id, None)
    except Exception as e:
        print(str(e), file=sys.stderr)
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

    # Chunk size needs to adjusted by how many genes are present, so that the payload always stays under the body size limit
    chunk_size = calculate_chunk_size(len(target_df.index), len(target_df.columns))

    print("TARGET: {}\nGENECART: {}\nTARGET DF (genes,samples): {}\nSAMPLES PER CHUNK: {}".format(dataset_id, genecart_id, target_df.shape, chunk_size), file=sys.stderr)

    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)

    try:
        results = loop.run_until_complete(fetch_all(loop, target_df, loading_df, is_pca, genecart_id, dataset_id, chunk_size))
    except Exception as e:
        print(str(e), file=sys.stderr)
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

    projection_patterns_df = concat_fetch_results_to_dataframe(results)

    # There is a good chance the samples are now out of order, which will break
    # the copying of the dataset observation metadata when this output is converted
    # to an AnnData object. So reorder back to dataset sample order.
    projection_patterns_df = projection_patterns_df.reindex(adata.obs.index.tolist())

    # Have had cases where the column names are x1, x2, x3, etc. so load in the original pattern names
    projection_patterns_df.set_axis(loading_df.columns, axis="columns", inplace=True)

    projection_patterns_df.to_csv(dataset_projection_csv)

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
    genecart_projection_csv = build_projection_csv_path(genecart_id, projection_id, "genecart")
    genecart_projection_json_file = build_projection_json_path(genecart_id, "genecart")

    # Symlink dataset_projection_csv to genecart_projection_csv (this syntax feels like it's in reverse)
    # NOTE: In the Docker instance, symlink reflects path to the mounted volume, not the local path
    try:
        genecart_projection_csv.symlink_to(dataset_projection_csv)
    except FileExistsError:
        print("Symlink already exists for {}".format(dataset_projection_csv), file=sys.stderr)

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

def run_tsne():
    if not gene_symbol or not dataset_id:
        return {
            "success": -1,
            "message": "Request needs both dataset id and gene symbol."
        }

    try:
        ana = get_analysis(analysis, dataset_id, session_id, analysis_owner_id)
    except PlotError as pe:
        return {
            "success": -1,
            "message": str(pe)
        }

    adata = ana.get_adata(backed=True)

    global projection_id
    if projection_id:
        try:
            adata = create_projection_adata(adata, dataset_id, projection_id)
        except PlotError as pe:
            return {
                'success': -1,
                'message': str(pe),
            }

    gene_symbols = (gene_symbol,)
    if 'gene_symbol' in adata.var.columns:
        gene_filter = adata.var.gene_symbol.isin(gene_symbols)
        if not gene_filter.any():
            return {
                'success': -1,
                'message': 'Gene not found',
            }

    else:
        return {
            'success': -1,
            'message': 'Missing gene_symbol in adata.var'
        }

    # Primary dataset - find tSNE_1 and tSNE_2 in obs and build X_tsne
    if analysis is None or analysis in ["null", "undefined", dataset_id]:
        for ds in [x_axis, y_axis]:
            if ds not in adata.obs:
                return {
                    'success': -1,
                    'message': 'Dataseries {} was selected but not present in adata.obs'.format(ds)
                }
        adata.obsm['X_tsne'] = adata.obs[[x_axis, y_axis]].values
        adata.obsm['X_umap'] = adata.obs[[x_axis, y_axis]].values
        adata.obsm['X_pca'] = adata.obs[[x_axis, y_axis]].values

    # We also need to change the adata's Raw var dataframe
    # We can't explicitly reset its index so we reinitialize it with
    # the newer adata object.
    # https://github.com/theislab/anndata/blob/master/anndata/base.py#L1020-L1022
    if adata.raw is not None:
        adata.raw = adata

    # Reorder the categorical values in the observation dataframe
    # Currently in UI only "plot_by_group" has reordering capabilities
    global order
    if order:
        order = order
        obs_keys = order.keys()
        for key in obs_keys:
            col = adata.obs[key]
            try:
                # Some columns might be numeric, therefore
                # we don't want to reorder these
                reordered_col = col.cat.reorder_categories(
                    order[key], ordered=True)
                adata.obs[key] = reordered_col
            except:
                pass

    selected = adata[:, gene_filter]
    df = selected.to_df()
    success = 1
    message = ""
    if len(df.columns) > 1:
        success = 2
        message = "WARNING: Multiple Ensemble IDs found for gene symbol '{}'.  Using the first stored Ensembl ID.".format(gene_symbol)

    # Drop duplicate gene symbols so that only 1 ensemble ID is used in scanpy
    adata.var = adata.var.reset_index().set_index('gene_symbol')
    # Currently the ensembl_id column is still called 'index', which could be confusing when looking at the new .index
    # Rename to end the confusion
    adata.var = adata.var.rename(columns={adata.var.columns[0]: "ensembl_id"})
    # Modify the AnnData object to not include any duplicated gene symbols (keep only first entry)
    if len(df.columns) > 1:
        scanpy_copy = ana.dataset_path().replace('.h5ad', '.scanpy_dups_removed.h5ad')
        if os.path.exists(scanpy_copy):
            os.remove(scanpy_copy)
        adata = adata[:, adata.var.index.duplicated() == False].copy(filename=scanpy_copy)

    io_fig = None
    try:
        basis = PLOT_TYPE_TO_BASIS[plot_type]
    except:
        raise("{} was not a valid plot type".format(plot_type))

    # NOTE: This may change in the future if users want plots by group w/o the colorize_by plot added
    if plot_by_group:
        global skip_gene_plot
        skip_gene_plot = None

    # Reverse cividis so "light" is at 0 and 'dark' is at incresing expression
    expression_color = create_colorscale_with_zero_gray("cividis_r" if colorblind_mode else "YlOrRd")

    # If colorize_by is passed we need to generate that image first, before the index is reset
    #  for gene symbols, then merge them.
    if colorize_by:
        # were custom colors passed?  the color index is the 'colorize_by' label but with '_colors' appended
        color_idx_name = "{0}_colors".format(colorize_by)

        ## why 2?  Handles the cases of a stringified "{}" or actual keyed JSON
        if colors is not None and len(colors) > 2:
            adata.uns[color_idx_name] = [colors[idx] for idx in adata.obs[colorize_by].cat.categories]

        elif color_idx_name in adata.obs:
            # Alternative method.  Associate with hexcodes already stored in the dataframe
            # Making the assumption that these values are hexcodes
            grouped = adata.obs.groupby([colorize_by, color_idx_name])
            # Ensure one-to-one mapping between category and hexcodes
            if len(adata.obs[colorize_by].unique()) == len(grouped):
                # Test if names are color hexcodes and use those if applicable (if first is good, assume all are)
                color_hex = adata.obs[color_idx_name].unique().tolist()
                if re.search(COLOR_HEX_PTRN, color_hex[0]):
                    color_map = {name[0]:name[1] for name, group in grouped}
                    adata.uns[color_idx_name] = [color_map[k] for k in adata.obs[colorize_by].cat.categories]


        if colorblind_mode:
            # build a cividis color map for the colorblind mode
            cb_colors = get_colorblind_scale(len(adata.obs[colorize_by].unique()))
            color_map = {name:cb_colors[idx] for idx, name in enumerate(adata.obs[colorize_by].cat.categories)}
            adata.uns[color_idx_name] = [color_map[k] for k in adata.obs[colorize_by].cat.categories]

        # Calculate the number of columns in the legend (if applicable)
        num_cols = calculate_num_legend_cols(len(adata.obs[colorize_by].unique()))

        # Get for legend order.
        colorize_by_order = adata.obs[colorize_by].unique()

        """
        NOTE: Quick note about legend "loc" and "bbox_to_anchor" attributes:

        bbox_to_anchor is the location of the legend relative to the plot frame.
        If x and y are 0, that is the lower-left corner of the plot.
        If bbox_to_anchor has 4 options, they are x, y, width, and height.  The last two are ratios relative to the plot. And x and y are the lower corner of the bounding box

        loc is the portion of the legend that will be at the bbox_to_anchor point.
        So, if x=0, y=0, and loc = "lower_left", the lower left corner of the legend will be anchored to the lower left corner of the plot
        """

        # If plotting by group the plot dimensions need to be determined
        if plot_by_group:
            column_order = adata.obs[plot_by_group].unique()
            group_len = len(column_order)
            num_plots = group_len + 2

            max_cols = num_plots
            if max_columns:
                max_cols = min(int(max_columns), num_plots)
            max_rows = ceil((num_plots) / max_cols)

            # Set up the figure specs
            figwidth = calculate_figure_width(max_cols)
            figheight = calculate_figure_height(max_rows)
            io_fig = plt.figure(figsize=(figwidth,figheight))
            spec = io_fig.add_gridspec(ncols=max_cols, nrows=max_rows)

            adata.obs["gene_expression"] = [float(x) for x in adata[:,adata.var.index.isin([gene_symbol])].X]
            max_expression = max(adata.obs["gene_expression"].tolist())

            row_counter = 0
            col_counter = 0

            # Filter expression data by "plot_by_group" group and plot each instance
            if order and plot_by_group in order:
                column_order = order[plot_by_group]

            for _,name in enumerate(column_order):
                # Copy gene expression dataseries to observation
                # Filter only expression values for a particular group.
                adata.obs["split_by_group"] = adata.obs.apply(lambda row: row["gene_expression"] if row[plot_by_group] == name else 0, axis=1)
                f = io_fig.add_subplot(spec[row_counter, col_counter])
                sc.pl.embedding(adata, basis=basis, color=["split_by_group"], color_map=expression_color, ax=f, show=False, use_raw=False, title=name, vmax=max_expression)
                col_counter += 1
                # Increment row_counter when the previous row is filled.
                if col_counter % max_cols == 0:
                    row_counter += 1
                    col_counter = 0
            # Add total gene plot and color plots
            if not skip_gene_plot:
                f_gene = io_fig.add_subplot(spec[row_counter, col_counter])    # final plot with colorize-by group
                sc.pl.embedding(adata, basis=basis, color=[gene_symbol], color_map=expression_color, ax=f_gene, show=False, use_raw=False) # Max expression is vmax by default
                col_counter += 1
                # Increment row_counter when the previous row is filled.
                if col_counter % max_cols == 0:
                    row_counter += 1
                    col_counter = 0
            f_color = io_fig.add_subplot(spec[row_counter, col_counter])    # final plot with colorize-by group
            sc.pl.embedding(adata, basis=basis, color=[colorize_by], ax=f_color, show=False, use_raw=False)
            (handles, labels) = sort_legend(f_color, colorize_by_order, horizontal_legend)
            f_color.legend(ncol=num_cols, bbox_to_anchor=[1, 1], frameon=False, handles=handles, labels=labels)
            if horizontal_legend:
                io_fig.legend(loc="upper center", bbox_to_anchor=[0, 0, 1, 0], frameon=False, ncol=NUM_HORIZONTAL_COLS, handles=handles, labels=labels)
                f_color.get_legend().remove()  # Remove legend added by scanpy
        else:
            # If 'skip_gene_plot' is set, only the colorize_by plot is printed, otherwise print gene symbol and colorize_by plots
            if skip_gene_plot:
                # the figsize options here (paired with dpi spec above) dramatically affect the definition of the image
                io_fig = plt.figure(figsize=(6, 4))
                if len(adata.obs[colorize_by].cat.categories) > 10:
                    io_fig = plt.figure(figsize=(13, 4))
                spec = io_fig.add_gridspec(ncols=1, nrows=1)
                f1 = io_fig.add_subplot(spec[0,0])
                sc.pl.embedding(adata, basis=basis, color=[colorize_by], ax=f1, show=False, use_raw=False)
                (handles, labels) = sort_legend(f1, colorize_by_order, horizontal_legend)
                f1.legend(ncol=num_cols, bbox_to_anchor=[1, 1], frameon=False, handles=handles, labels=labels)
                if horizontal_legend:
                    io_fig.legend(loc="upper center", bbox_to_anchor=[0, 0, 1, 0], frameon=False, ncol=NUM_HORIZONTAL_COLS, handles=handles, labels=labels)
                    f1.get_legend().remove()  # Remove legend added by scanpy
            else:
                # the figsize options here (paired with dpi spec above) dramatically affect the definition of the image
                io_fig = plt.figure(figsize=(13, 4))
                spec = io_fig.add_gridspec(ncols=2, nrows=1, width_ratios=[1.1, 1])
                f1 = io_fig.add_subplot(spec[0,0])
                f2 = io_fig.add_subplot(spec[0,1])
                sc.pl.embedding(adata, basis=basis, color=[gene_symbol], color_map=expression_color, ax=f1, show=False, use_raw=False)
                sc.pl.embedding(adata, basis=basis, color=[colorize_by], ax=f2, show=False, use_raw=False)
                (handles, labels) = sort_legend(f2, colorize_by_order, horizontal_legend)
                f2.legend(ncol=num_cols, bbox_to_anchor=[1, 1], frameon=False, handles=handles, labels=labels)
                if horizontal_legend:
                    io_fig.legend(loc="upper center", bbox_to_anchor=[0, 0, 1, 0], frameon=False, ncol=NUM_HORIZONTAL_COLS, handles=handles, labels=labels)
                    f2.get_legend().remove()  # Remove legend added by scanpy
    else:
        io_fig = sc.pl.embedding(adata, basis=basis, color=[gene_symbol], color_map=expression_color, return_fig=True, use_raw=False)

    io_pic = io.BytesIO()
    io_fig.tight_layout()   # This crops out much of the whitespace around the plot. The next line does this with the legend too
    io_fig.savefig(io_pic, format='png', bbox_inches="tight")
    io_pic.seek(0)
    plt.clf()
    plt.close()  # Prevent zombie plots, which can cause issues

    return {
        "success": success,
        "message": message,
        "image": base64.b64encode(io_pic.read()).decode("utf-8")
    }

def main():
    print("running prestep", file=sys.stderr)
    run_projection_prestep()
    print("running projection", file=sys.stderr)
    res = run_projection()
    if res.get("success", 0) == 0:
        print(res.get("message", "No error message"))
        sys.exit(1)

if __name__ == "__main__":
    main()
    print("running tSNE", file=sys.stderr)
    res2 = run_tsne()
    if res2.get("success", 0) == 0:
        print(res2.get("message", "No error message"))
        sys.exit(1)
    print("Successful run", file=sys.stderr)
    sys.exit(0)