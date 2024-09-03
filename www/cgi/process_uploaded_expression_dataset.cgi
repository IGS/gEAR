#!/opt/bin/python3

"""
Process the uploaded expression dataset, regardless of type.  As the data are,
a process which can take hours, the following data structure is periodically 
written to the same directory as the dataset:

status.json

{
    "process_id": 1234,
    "status": "processing",
    "message": "Processing the dataset.  This may take a while.",
    "progress": 0
}

Where status can be 'extracting', 'processing', 'error', or 'complete'.
"""

import cgi
import json
import os, sys

import pandas as pd
import scanpy as sc
from scipy import sparse

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
    "process_id": os.getpid(),
    "status": "extracting",
    "message": "",
    "progress": 0
}

def main():
    result = {'success':0 }
    global share_uid
    global session_id

    form = cgi.FieldStorage()
    share_uid = form.getvalue('share_uid')
    session_id = form.getvalue('session_id')

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        print_and_go(None, json.dumps(result))

    # values are mex_3tab, excel, rdata, h5ad
    dataset_formats = ['mex_3tab', 'excel', 'rdata', 'h5ad']
    dataset_format = form.getvalue('dataset_format')
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

    # Since this process can take a while, we want to fork off of apache and continue
    #  processing in the background.
    with open(status_file, 'w') as f:
        f.write(json.dumps(status))

    ###############################################
    # This is the fork off of apache
    # https://stackoverflow.com/a/22181041/1368079
    sys.stdout = original_stdout
    result['success'] = 1
    print(json.dumps(result)
    
    sys.stdout.flush()
    os.close(sys.stdout.fileno()) # Break web pipe
    sys.stderr.flush()
    os.close(sys.stderr.fileno()) # Break web pipe

    if os.fork(): # Get out parent process
        if dataset_format == 'mex_3tab':
            process_mex_3tab(dataset_upload_dir)

    sys.exit(0)


def print_and_go(status_file, content):
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

    status['status'] = 'processing'
    status['message'] = 'Initializing dataset processing.'
    with open(os.path.join(upload_dir, 'status.json'), 'w') as f:
        f.write(json.dumps(status))

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
    for chunk in reader:
        rows_read += chunk_size
        percentage = int((rows_read / total_rows) * 100)
        expression_matrix.append(sparse.csr_matrix(chunk.values))
        
        status['progress'] = percentage
        status['message'] = f"Processed {rows_read}/{total_rows} expression matrix chunks ..."
        with open(os.path.join(upload_dir, 'status.json'), 'w') as f:
            f.write(json.dumps(status))

    adata.X = sparse.vstack(expression_matrix)
    
    print("Finished reading expression matrix file", file=sys.stderr, flush=True)
    adata = adata.transpose()

    h5ad_path = os.path.join(upload_dir, f"{share_uid}.h5ad")
    adata.write(h5ad_path)

    status['status'] = 'complete'
    status['message'] = 'Dataset processed successfully.'
    with open(os.path.join(upload_dir, 'status.json'), 'w') as f:
        f.write(json.dumps(status))
    

def process_mex(upload_dir):
    pass

def process_mex_3tab(upload_dir):
    # Extract the file
    import tarfile
    filename = os.path.join(upload_dir, f"{share_uid}.tar.gz")

    files_extracted = []

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

    # Determine the dataset type
    dataset_type = tarball_content_type(files_extracted)

    if dataset_type is None:
        status['status'] = 'error'
        status['message'] = "Unsupported dataset format. Couldn't tell type from file names within the tarball"
        with open(os.path.join(upload_dir, 'status.json'), 'w') as f:
            f.write(json.dumps(status))

    # Call the appropriate function
    if dataset_type == 'threetab':
        process_3tab(upload_dir)
    elif dataset_type == 'mex':
        process_mex(upload_dir)

def tarball_content_type(filenames):
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
