#!/opt/bin/python3

import cgi
import json
import os, sys

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDERR, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

from gear.datasetuploader import FileType, DatasetUploader

# This is another attempt to fix the STDOUT hack above in a cleaner way, but I couldn't
#  get it to actually change the loglevel of another imported module.
#import logging
#logging.getLogger("matplotlib").setLevel(logging.WARNING)

def main():
    user_upload_file_base = '../uploads/files'

    result = {'success':0, 'shape': {}, 'preview_gene': '', 'plots': [] }

    form = cgi.FieldStorage()
    filename = form.getvalue('filename')

    # This means the PHP upload failed
    if filename is None:
        result['message'] = "File upload failed.  Try again and contact us if this continues."
        print_and_go(json.dumps(result))

    filepath = user_upload_file_base + "/" + filename
    h5_out_path = user_upload_file_base + "/" + filename.replace('.', '_') + '.h5ad'
    session_id = form.getvalue('session_id')

    user = geardb.get_user_from_session_id(session_id)

    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
    else:
        # gEAR User found. Read and Validate expression file
        file_type = filename.rsplit('.', 1)[1]

        #Initize uploader for expression data (with factory)
        dsu = DatasetUploader()
        dataset_uploader = dsu.get_by_filetype(filetype=file_type, filepath=filepath)
        try:
            dataset_uploader._read_file(filepath)

            #Let's give the user some feedback on the shape of the dataset.
            rows_X, cols_X = dataset_uploader.adata.X.shape
            cols_obs, rows_obs = dataset_uploader.adata.obs.shape
            cols_var, rows_var = dataset_uploader.adata.var.shape

            dataset_uploader.adata.write(h5_out_path)

            result['success'] = 1
            result['message'] = 'File successfully read.'

            result['shape']['X'] = {'rows': str(rows_X), 'cols': str(cols_X)}
            result['shape']['obs'] = {'rows': str(rows_obs), 'cols': str(cols_obs)}
            result['shape']['var'] = {'rows': str(rows_var), 'cols': str(cols_var)}

        except Exception as err:
            result['message'] = str(err)

    print_and_go(json.dumps(result))

def print_and_go(content):
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n', flush=True)
    print(content)
    sys.exit(0)

if __name__ == '__main__':
    main()
