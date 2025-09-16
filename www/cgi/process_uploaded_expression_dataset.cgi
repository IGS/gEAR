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
import time
import typing
import zipfile

import anndata
import pandas as pd
import scanpy as sc
from scipy import sparse

if typing.TYPE_CHECKING:
    from typing import NoReturn



# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDOUT, killing this CGI.  So here we redirect
#  STDOUT until we need it.
print('Content-Type: application/json\n\n', flush=True)
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

share_uid = None
session_id = None
user_upload_file_base = '../uploads/files'

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


    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        print_and_go(None, json.dumps(result))

    # values are mex_3tab, excel, rdata, h5ad
    dataset_formats = ['mex_3tab', 'excel', 'rdata', 'h5ad', 'spatial']
    dataset_upload_dir = os.path.join(user_upload_file_base, session_id, share_uid)

    # quickly write the status so the page doesn't error out
    status_file = os.path.join(dataset_upload_dir, 'status.json')
    with open(status_file, 'w') as f:
        f.write(json.dumps(status))

    # if the upload directory doesn't exist, we can't process the dataset
    if not os.path.exists(dataset_upload_dir):
        result['message'] = 'Dataset/directory not found.'
        print_and_go(status_file, json.dumps(result))

    if dataset_format not in dataset_formats:
        result['message'] = 'Unsupported dataset format.'
        print_and_go(status_file, json.dumps(result))

    if dataset_format == "spatial":
        from gear.spatialhandler import SPATIALTYPE2CLASS
        if spatial_format not in SPATIALTYPE2CLASS:
            result['message'] = 'Invalid spatial format specified.'
            print_and_go(status_file, json.dumps(result))

    # Since this process can take a while, we want to fork off of apache and continue
    #  processing in the background.
    with open(status_file, 'w') as f:
        f.write(json.dumps(status))

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
    elif dataset_format == "spatial":
        spatial_handler_class = SPATIALTYPE2CLASS[spatial_format]
        """
            sp_class._read_file(args.input_file, organism_id=args.organism_id, dataset_id=args.dataset_id)
            output_filename = args.dataset_id
            output_path = DEST_DIRPATH / (output_filename + OUTPUT_SUFFIX)
            print("Writing to {0}".format(output_path))
            sp_class.write_to_zarr(filepath=output_path)
        """

        spatial_handler = spatial_handler_class(share_uid, dataset_upload_dir, status_file)
        spatial_handler.process_spatial()
    else:
        raise Exception('Unsupported dataset format')


def print_and_go(status_file, content) -> "NoReturn":
    sys.stdout = original_stdout
    print(content)

    if status_file is not None:
        with open(status_file, 'w') as f:
            f.write(json.dumps(status))

    sys.exit(0)

def process_3tab(upload_dir):
    import subprocess

    chunk_size = 500
    expression_matrix_path = None
    obs = None
    var = None

    write_status(upload_dir, 'processing', 'Initializing dataset processing.')

    for infile in os.listdir(upload_dir):
        # skip any files beginning with a dot
        if infile.startswith('.'):
            continue

        filepath = "{0}/{1}".format(upload_dir, infile)

        # Read each file as pandas dataframes
        if infile == 'expression.tab' or os.path.basename(filepath) == 'expression.tab' or 'DataMTX.tab' in infile:
            expression_matrix_path = filepath
        elif infile == 'observations.tab' or os.path.basename(filepath) == 'observations.tab' or 'COLmeta.tab' in infile:
            #print("Reading observations file: {0}".format(filepath), file=sys.stderr, flush=True)
            obs = pd.read_table(filepath, sep='\t', index_col=0, header=0)
        elif infile == 'genes.tab' or os.path.basename(filepath) == 'genes.tab' or 'ROWmeta.tab' in infile:
            #print("Reading genes file: {0}".format(filepath), file=sys.stderr, flush=True)
            var = pd.read_table(filepath, sep='\t', index_col=0, header=0)

    for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
        if str_type in obs.columns:
            obs[str_type] = pd.Categorical(obs[str_type])

    for num_type in ['replicate', 'time_point_order']:
        if num_type in obs.columns:
            obs[num_type] = pd.to_numeric(obs[num_type])

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
            with open(os.path.join(upload_dir, 'status.json'), 'w') as f:
                f.write(json.dumps(status))

        adata.X = sparse.vstack(expression_matrix)
    except Exception as e:
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
                with open(os.path.join(upload_dir, 'status.json'), 'w') as f:
                    f.write(json.dumps(status))

            except Exception as inner_e:
                #print(f"\nError in chunk {chunk_index}: {inner_e}")
                #print("Chunk head:")
                #print(chunk.head())
                raise

        # Try stacking the cleaned chunks
        try:
            adata.X = sparse.vstack(expression_matrix)
        except Exception as final_e:
            #print(f"\nFinal vstack still failed: {final_e}")

            #print("Collected chunk shapes:")
            #for i, shape in enumerate(chunk_shapes):
            #    print(f"  Chunk {i+1}: {shape}")

            raise


    adata = adata.transpose()
    adata.obs = sanitize_obs_for_h5ad(adata.obs)

    h5ad_path = os.path.join(upload_dir, f"{share_uid}.h5ad")
    adata.write(h5ad_path)

    write_status(upload_dir, 'complete', 'Dataset processed successfully.')

def process_excel(upload_dir):
    filepath = os.path.join(upload_dir, f"{share_uid}.xlsx")

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

    for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
        if str_type in obs_df.columns:
            obs_df[str_type] = pd.Categorical(obs_df[str_type])

    for num_type in ['replicate', 'time_point_order']:
        if num_type in obs_df.columns:
            obs_df[num_type] = pd.to_numeric(obs_df[num_type])

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

    h5ad_path = os.path.join(upload_dir, f"{share_uid}.h5ad")
    adata.write(h5ad_path)

    write_status(upload_dir, 'complete', 'Dataset processed successfully.')

def process_mex(upload_dir):
    pass

def process_mex_3tab(upload_dir):
    # Extract the file
    import tarfile
    compression_format = None
    filename = os.path.join(upload_dir, f"{share_uid}.tar.gz")

    if os.path.exists(filename):
        compression_format = 'tarball'
    else:
        filename = os.path.join(upload_dir, f"{share_uid}.zip")

        if os.path.exists(filename):
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
                            os.rename(os.path.join(upload_dir, entry.name),
                                    os.path.join(upload_dir, suffix))

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
                            os.rename(os.path.join(upload_dir, entry.filename),
                                    os.path.join(upload_dir, suffix))

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

def sanitize_obs_for_h5ad(obs_df):
    for col in obs_df.columns:
        if obs_df[col].dtype == 'object':
            obs_df[col] = obs_df[col].fillna('').astype(str)
    return obs_df

def write_status(upload_dir, status_name, message):
    status['status'] = status_name
    status['message'] = message
    with open(os.path.join(upload_dir, 'status.json'), 'w') as f:
        f.write(json.dumps(status))

def package_content_type(filenames):
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


if __name__ == '__main__':
    main()
