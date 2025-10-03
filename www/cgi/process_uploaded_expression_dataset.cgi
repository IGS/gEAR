#!/opt/bin/python3

"""
Process the uploaded expression dataset, regardless of type.  As the data are uploaded,
a process which can take hours, the following data structure is periodically
written to the same directory as the dataset:

status.json

{
    "process_id": 1234,
    "status": "processing",
    "message": "Processing the dataset.  This may take a while.",
    "progress": 0
}

Where status can be 'uploaded', 'extracting', 'processing', 'error', or 'complete'.
"""

import cgi
import json
import os
import sys
import zipfile
from pathlib import Path

import anndata
import pandas as pd
import scanpy as sc
from scipy import sparse

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDOUT, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = Path(__file__).resolve().parents[2] / 'lib'
sys.path.append(str(lib_path))
import geardb
from gear.spatialhandler import SPATIALTYPE2CLASS

share_uid = None
session_id = None
user_upload_file_base = '../uploads/files'

# Initial status (will be updated as we progress)
status = {
    "process_id": None,
    "status": "extracting",
    "message": "",
    "progress": 0
}

def main():
    result = {'success':0, "message": ""}
    global share_uid
    global session_id

    form = cgi.FieldStorage()
    share_uid = form.getvalue('share_uid')
    session_id = form.getvalue('session_id')
    dataset_format = form.getvalue('dataset_format')
    spatial_format = form.getvalue('spatial_format')  # may be None

    if share_uid is None or session_id is None or dataset_format is None:
        result['message'] = 'Missing one or more required parameters.'
        return result

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        return result

    # values are mex_3tab, excel, rdata, h5ad
    dataset_formats = ['mex_3tab', 'excel', 'rdata', 'h5ad', 'spatial']
    dataset_upload_dir = Path(user_upload_file_base) / session_id / share_uid

    # quickly write the status so the page doesn't error out

    status_file = Path(dataset_upload_dir) / 'status.json'
    with open(status_file, 'w') as f:
        f.write(json.dumps(status))

    # if the upload directory doesn't exist, we can't process the dataset
    if not dataset_upload_dir.is_dir():
        result['message'] = 'Dataset/directory not found.'
        write_status(dataset_upload_dir, 'error', result['message'])
        return result

    if dataset_format not in dataset_formats:
        result['message'] = 'Unsupported dataset format.'
        write_status(dataset_upload_dir, 'error', result['message'])
        return result

    if dataset_format == "spatial":
        if spatial_format not in SPATIALTYPE2CLASS:
            result['message'] = 'Invalid spatial format specified.'
            write_status(dataset_upload_dir, 'error', result['message'])
            return result

    # Since this process can take a while, we want to fork off of apache and continue
    #  processing in the background.

    ###############################################
    # This is the fork off of apache
    # https://stackoverflow.com/a/22181041/1368079
    # https://stackoverflow.com/questions/6024472/start-background-process-daemon-from-cgi-script
    # https://groups.google.com/g/comp.lang.python/c/gSRnd0RoVKY?pli=1
    do_fork = False
    if do_fork:
        sys.stdout = original_stdout
        result['success'] = 1
        print(json.dumps(result))

        sys.stdout.flush()
        sys.stdout.close()

        sys.stderr.flush()
        sys.stderr.close()

        f = os.fork()

        if f != 0:
            # Terminate the parent
            sys.exit(0)
            # PARENT DEAD

        # CHILD CONTINUES FROM HERE

    status['process_id'] = os.getpid()

    # new child command
    if dataset_format == 'mex_3tab':
        process_mex_3tab(dataset_upload_dir)
    elif dataset_format == 'excel':
        process_excel(dataset_upload_dir)
    elif dataset_format == "h5ad":
        process_h5ad(dataset_upload_dir)
    elif dataset_format == "spatial":
        process_spatial(dataset_upload_dir, spatial_format)
    else:
        result["success"] = 0
        result["message"] = f"Unsupported dataset format: {dataset_format}"
    return result

def process_h5ad(upload_dir: Path) -> None:
    """
    Processes an uploaded .h5ad (AnnData) file in the specified upload directory by performing the following steps:
    1. Reads the .h5ad file as an AnnData object.
    2. Sanitizes and categorizes the observation (obs) dataframe columns.
    3. Writes the sanitized data to a new .h5ad file.
    4. Replaces the original file with the sanitized version.
    5. Updates the processing status at each stage.

    Args:
        upload_dir (Path): The directory containing the uploaded .h5ad file.

    Returns:
        None
    """
    # If the file is an h5ad, it should be formatted as an AnnData object already.
    # But we still want to do some sanitization of the obs dataframe.

    write_status(upload_dir, 'processing', 'Initializing dataset processing.')

    filepath = upload_dir / f"{share_uid}.h5ad"
    adata = anndata.read_h5ad(filepath)
    obs = adata.obs

    categorize_observation_columns(obs)
    adata.obs = sanitize_obs_for_h5ad(obs)

    h5ad_path = upload_dir / f"{share_uid}.new.h5ad"
    adata.write(h5ad_path)

    # Replace the original file with the sanitized one
    filepath.unlink()  # remove original
    h5ad_path.rename(filepath)  # rename new to original name

    write_status(upload_dir, 'complete', 'Dataset processed successfully.')

def process_3tab(upload_dir: Path) -> None:
    import subprocess

    chunk_size = 500
    expression_matrix_path = None
    obs = None
    var = None

    write_status(upload_dir, 'processing', 'Initializing dataset processing.')

    for infile in upload_dir.iterdir():
        # skip any files beginning with a dot
        if infile.name.startswith('.'):
            continue

        filepath = "{0}/{1}".format(upload_dir, infile)

        # Read each file as pandas dataframes
        if infile.name == 'expression.tab' or 'DataMTX.tab' in infile.name:
            expression_matrix_path = filepath
        elif infile.name == 'observations.tab' or 'COLmeta.tab' in infile.name:
            #print("Reading observations file: {0}".format(filepath), file=sys.stderr, flush=True)
            obs = pd.read_table(filepath, sep='\t', index_col=0, header=0)
        elif infile.name == 'genes.tab' or 'ROWmeta.tab' in infile.name:
            #print("Reading genes file: {0}".format(filepath), file=sys.stderr, flush=True)
            var = pd.read_table(filepath, sep='\t', index_col=0, header=0)

    if obs is None:
        write_status(upload_dir, 'error', "No observations file found. Expected observations.tab or COLmeta.tab.")
        return

    if var is None:
        write_status(upload_dir, 'error', "No genes file found. Expected genes.tab or ROWmeta.tab.")
        return

    if expression_matrix_path is None:
        write_status(upload_dir, 'error', "No expression file found. Expected expression.tab or DataMTX.tab.")
        return

    categorize_observation_columns(obs)

    # Read in expressions as AnnData object in a memory-efficient manner
    adata = sc.AnnData(obs=var, var=obs)
    reader = pd.read_csv(expression_matrix_path, sep='\t', index_col=0, chunksize=chunk_size)

    # This can be an order of magnitude faster than the using python alone
    total_rows = int(subprocess.check_output(f"/usr/bin/wc -l {expression_matrix_path}", shell=True).split()[0])

    expression_matrix = []
    rows_read = 0

    ## Try to process the file the quickest way first, assuming things are peachy. Then, if not,
    #  do some checks and conversions (slower) as a backup.
    try:
        for chunk in reader:
            rows_read += chunk_size
            percentage = int((rows_read / total_rows) * 100)
            expression_matrix.append(sparse.csr_matrix(chunk.values))

            status['progress'] = percentage
            status['message'] = f"Processed {rows_read}/{total_rows} expression matrix chunks ..."
            with open(upload_dir / "status.json", 'w') as f:
                f.write(json.dumps(status))

        adata.X = sparse.vstack(expression_matrix) # type: ignore
    except Exception:
        #print(f"\nOriginal vstack failed: {e}")
        #print("Retrying with per-chunk cleanup...")

        expression_matrix.clear()
        chunk_shapes = []
        rows_read = 0
        status['progress'] = 0

        # Re-open reader here
        reader = pd.read_csv(expression_matrix_path, sep='\t', index_col=0, chunksize=chunk_size)

        for chunk_index, chunk in enumerate(reader, start=1):
            try:
                 # Clean each cell: strip string values
                chunk = chunk.replace(r'^\s+|\s+$', '', regex=True)
                chunk = chunk.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))

                # Convert to numeric (non-numeric → NaN → fill with 0)
                chunk_numeric = chunk.apply(pd.to_numeric, errors='coerce').fillna(0)

                matrix = sparse.csr_matrix(chunk_numeric.values)
                expression_matrix.append(matrix)
                chunk_shapes.append(matrix.shape)

                rows_read += chunk_size
                percentage = int((rows_read / total_rows) * 100)

                status['progress'] = percentage
                status['message'] = f"Processed {rows_read}/{total_rows} expression matrix chunks ..."
                with open(upload_dir / "status.json", 'w') as f:
                    f.write(json.dumps(status))

            except Exception:
                #print(f"\nError in chunk {chunk_index}: {inner_e}")
                #print("Chunk head:")
                #print(chunk.head())
                raise

        # Try stacking the cleaned chunks
        try:
            adata.X = sparse.vstack(expression_matrix) # type: ignore
        except Exception:
            #print(f"\nFinal vstack still failed: {final_e}")

            #print("Collected chunk shapes:")
            #for i, shape in enumerate(chunk_shapes):
            #    print(f"  Chunk {i+1}: {shape}")

            raise


    adata = adata.transpose()
    adata.obs = sanitize_obs_for_h5ad(adata.obs)

    h5ad_path = upload_dir / f"{share_uid}.h5ad"
    adata.write(h5ad_path)

    write_status(upload_dir, 'complete', 'Dataset processed successfully.')

def process_excel(upload_dir: Path) -> None:
    filepath = upload_dir / f"{share_uid}.xlsx"

    write_status(upload_dir, 'processing', 'Initializing dataset processing.')

    exp_df = pd.read_excel(filepath, sheet_name='expression', index_col=0).transpose()

    try:
        X = exp_df.to_numpy()[:, 0:].astype(float)
    except ValueError:
        write_status(upload_dir, 'error', "Encountered unexpected value type. Expected float type in expression matrix.")
        return

    # Get counts of genes and observations
    number_obs_from_exp, number_genes_from_exp = X.shape

    # Get the observations
    try:
        obs_df = pd.read_excel(filepath, sheet_name='observations', index_col='observations')
    except ValueError:
        write_status(upload_dir, 'error', "No observations sheet found. Expected spreadsheet sheet named 'observations'.")
        return

    # Verify number observations equal those found in expression sheet
    number_obs, number_cond = obs_df.shape
    if number_obs != number_obs_from_exp:
        write_status(upload_dir, 'error', "Observations sheet error. Row count ({0}) in 'observations' sheet must match row count of 'expression' sheet({1}).".format(number_obs, number_obs_from_exp))
        return

    # Verify observations index matches expression sheet col index
    if not obs_df.index.equals(exp_df.index):
        write_status(upload_dir, 'error', "Observations sheet error. The names and order of the index column in 'observations' sheet must match the rows of 'expression' sheet.")
        return

    # Get the genes (if present), else use the .var from exp_df
    try:
        genes_df = pd.read_excel(filepath, sheet_name='genes', index_col=0, converters={'gene_symbol': str})
    except ValueError:
        # With 'genes' sheet absent try to get the genes from 'expression'
        try:
            # read expression sheet. Set genes as the index and only parse 1 column
            genes_df = pd.read_excel(filepath, sheet_name='expression', index_col=0, usecols=[0,1])

            # remove the 1 parsed column. This leaves an empty dataframe
            # with an index of gene ids to match other datasets.
            genes_df = genes_df.drop(genes_df.columns[0], axis=1)
        except Exception as err:
            write_status(upload_dir, 'error', "No 'genes' sheet found. Tried using genes column from 'expression' sheet as .var, but " + str(err))
            return

    # Check for numeric gene symbols
    if 'gene_symbol' in genes_df.columns:
        digit_count = genes_df['gene_symbol'].str.isnumeric().sum()
        if digit_count > 0:
            write_status(upload_dir, 'error', "Genes sheet error. {0} gene symbols are listed as numbers, not gene symbols".format(str(digit_count)))
            return
    else:
        write_status(upload_dir, 'error', "Failed to find gene_symbol column in genes tab")
        return

    categorize_observation_columns(obs_df)

    # Verify number genes equal those found in expression sheet
    number_genes, number_columns = genes_df.shape
    if number_genes != number_genes_from_exp:
        write_status(upload_dir, 'error', "Genes sheet error. Row count ({0}) in 'genes' sheet must match row count of 'expression' sheet({1}).".format(number_genes, number_genes_from_exp))
        return

    # Verify genes index matches expression sheet columns
    if not genes_df.index.equals(exp_df.columns):
        write_status(upload_dir, 'error', "Genes sheet error. The names and order of the index column in 'genes' sheet must match the rows of 'expression' sheet.")
        return

    # Create AnnData object and return it
    adata = anndata.AnnData(X=X, obs=obs_df, var=genes_df)
    adata.obs = sanitize_obs_for_h5ad(adata.obs)

    h5ad_path = upload_dir / f"{share_uid}.h5ad"
    adata.write(h5ad_path)

    write_status(upload_dir, 'complete', 'Dataset processed successfully.')

def process_mex(upload_dir: Path) -> None:
    pass

def process_mex_3tab(upload_dir: Path) -> None:
    # Extract the file
    import tarfile
    compression_format = None
    filename = upload_dir / f"{share_uid}.tar.gz"

    if filename.exists():
        compression_format = 'tarball'
    else:
        filename = upload_dir / f"{share_uid}.zip"

        if filename.exists():
            compression_format = 'zip'
        else:
            write_status(upload_dir, 'error', "No tarball or zip file found.")
            return

    files_extracted = []

    if compression_format == 'tarball':
        try:
            with tarfile.open(filename) as tf:
                for entry in tf:
                    tf.extract(entry, path=upload_dir)

                    # Nemo suffixes
                    nemo_suffixes = ['DataMTX.tab', 'COLmeta.tab', 'ROWmeta.tab']
                    suffix_found = None

                    for suffix in nemo_suffixes:
                        if entry.name.endswith(suffix):
                            suffix_found = suffix
                            # Rename the file to the appropriate name
                            old_name = upload_dir / entry.name
                            new_name = upload_dir / suffix
                            old_name.replace(new_name)

                    if suffix_found is not None:
                        files_extracted.append(suffix_found)
                    else:
                        files_extracted.append(entry.name)
        except tarfile.ReadError:
            write_status(upload_dir, 'error', "Bad tarball file. Couldn't extract the tarball.")
            return

    if compression_format == 'zip':
        try:
            with zipfile.ZipFile(filename) as zf:
                for entry in zf.infolist():
                    zf.extract(entry, path=upload_dir)

                    # Nemo suffixes
                    nemo_suffixes = ['DataMTX.tab', 'COLmeta.tab', 'ROWmeta.tab']
                    suffix_found = None

                    for suffix in nemo_suffixes:
                        if entry.filename.endswith(suffix):
                            suffix_found = suffix
                            # Rename the file to the appropriate name
                            old_name = upload_dir / entry.filename
                            new_name = upload_dir / suffix
                            old_name.replace(new_name)

                    if suffix_found is not None:
                        files_extracted.append(suffix_found)
                    else:
                        files_extracted.append(entry.filename)
        except zipfile.BadZipFile:
            write_status(upload_dir, 'error', "Bad zip file. Couldn't extract the zip file.")
            return

    # Determine the dataset type
    dataset_type = package_content_type(files_extracted)

    if dataset_type is None:
        write_status(upload_dir, 'error', "Unsupported dataset format. Couldn't tell type from file names within the tarball")

    # Call the appropriate function
    if dataset_type == 'threetab':
        process_3tab(upload_dir)
    elif dataset_type == 'mex':
        process_mex(upload_dir)

def process_spatial(upload_dir: Path, spatial_format: str) -> None:
    """
    Processes a spatial transcriptomics dataset uploaded to a specified directory.

    This function handles the reading and conversion of spatial data files using a handler
    class determined by the spatial_format. It expects a metadata.json file in the upload
    directory to extract the sample's taxonomic ID, which is then used to retrieve the
    organism ID. The function reads the spatial data archive, processes it, and writes the
    output in Zarr format. Status updates and errors are logged using the write_status function.

    Args:
        upload_dir (str): The directory where the uploaded files are located.
        spatial_format (str): The format of the spatial data, used to select the appropriate handler.

    Raises:
        Writes error status if the metadata file is missing or if reading/converting the spatial file fails.
    """
    spatial_obj = SPATIALTYPE2CLASS[spatial_format]()   # instantiate the appropriate handler class
    metadata_file = upload_dir / 'metadata.json'
    if not metadata_file.is_file():
        write_status(upload_dir, 'error', "No metadata JSON file found.")

    # get organism_id by converting sample_taxid(needed for some but not all spatial handlers)
    with open(metadata_file, 'r') as f:
        metadata = json.load(f)
    sample_taxid = metadata.get("sample_taxid", None)
    organism_id=geardb.get_organism_id_by_taxon_id(sample_taxid)
    filepath = upload_dir / f"{share_uid}.tar.gz"

    try:
        spatial_obj.process_file(filepath, extract_dir=upload_dir, organism_id=organism_id)
    except Exception as e:
        write_status(upload_dir, 'error', f"Error in uploading spatial file: {e}")
        return

    output_filename = f"{share_uid}.zarr"
    output_path = upload_dir / output_filename
    # Remove existing Zarr store if it exists
    # This is a safeguard; it shouldn't exist at this point, unless there was a failure post-writing
    if output_path.exists():
        import shutil
        shutil.rmtree(output_path)

    write_status(upload_dir, 'processing', 'Writing Zarr store')
    spatial_obj.write_to_zarr(filepath=output_path)
    write_status(upload_dir, 'complete', 'Dataset processed successfully.')

def sanitize_obs_for_h5ad(obs_df: pd.DataFrame) -> pd.DataFrame:
    for col in obs_df.columns:
        if obs_df[col].dtype == 'object':
            obs_df[col] = obs_df[col].fillna('').astype(str)
    return obs_df

def write_status(upload_dir: Path, status_name: str, message: str) -> None:
    status['status'] = status_name
    status['message'] = message
    with open(upload_dir / 'status.json', 'w') as f:
        f.write(json.dumps(status))

def package_content_type(filenames: list[str]) -> str | None:
        #print("DEBUG: filenames", file=sys.stderr, flush=True)
        #print(filenames, file=sys.stderr, flush=True)
        """
        mex:
        matrix.mtx
        barcodes.tsv
        genes.tsv

        threetab:
        expression.tab
        genes.tab
        observations.tab

        None is returned if neither of these is true

        Added NEMO file format functionality.
        DataMTX.tab -> expression.tab
        COLmeta.tab -> observations.tab
        ROWmeta.tab -> genes.tab
        """
        if 'expression.tab' in filenames and 'genes.tab' in filenames and 'observations.tab' in filenames:
            return 'threetab'

        if 'matrix.mtx' in filenames and 'barcodes.tsv' in filenames and 'genes.tsv' in filenames:
            return 'mex'

        if 'DataMTX.tab' in filenames and 'COLmeta.tab' in filenames and 'ROWmeta.tab' in filenames:
            return 'threetab'

        return None

def categorize_observation_columns(obs: pd.DataFrame) -> None:
    """
    Categorizes and converts specific columns in the observation DataFrame.

    For each of the string-type columns ('cell_type', 'condition', 'time_point', 'time_unit') present in the DataFrame,
    converts the column to a pandas Categorical type. For each of the numeric-type columns ('replicate', 'time_point_order')
    present in the DataFrame, converts the column to a numeric type.

    Parameters:
        obs (pd.DataFrame): The observation DataFrame whose columns will be categorized or converted.

    Returns:
        None: The function modifies the DataFrame in place.
    """
    for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
        if str_type in obs.columns:
            obs[str_type] = pd.Categorical(obs[str_type])

    for num_type in ['replicate', 'time_point_order']:
        if num_type in obs.columns:
            obs[num_type] = pd.to_numeric(obs[num_type])

if __name__ == '__main__':
    result = main()
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n', flush=True)
    print(json.dumps(result), flush=True)
