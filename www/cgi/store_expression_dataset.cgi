#!/opt/bin/python3

"""
Used by the expression uploader, this stores the actual dataset from the form
and saves it to a file for processing.

Writes a file at: ../uploads/files/<session_id>/<share_uid>/<share_uid>.<ext>
"""

import cgi
import json
import sys
from pathlib import Path

lib_path = Path(__file__).resolve().parents[2] / 'lib'
sys.path.append(str(lib_path))
import geardb

def main():
    print('Content-Type: application/json\n\n')
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    share_uid = form.getvalue('share_uid')
    dataset_format = form.getvalue('dataset_format')
    spatial_format = form.getvalue('spatial_format')  # may be None

    if not share_uid: # should never happen
        error_msg = f"Unexpected missing share_uid in store_expression_dataset.cgi. session_id={session_id!r}"
        print(error_msg, file=sys.stderr)
        result = {'success': 0, 'message': 'Internal error: share_uid missing (this should never happen). Please contact support.'}
        return result

    user = geardb.get_user_from_session_id(session_id)
    result = {'success': 0, 'message': ''}

    filename = form['dataset_file'].filename

    if filename.endswith('.tar.gz'):
        file_extension = 'tar.gz'
    else:
        file_extension = filename.split('.')[-1]

    # This should already have been created when the metadata was stored
    user_upload_file_base = "../uploads/files/{0}".format(session_id)

    dataset_filename = Path(user_upload_file_base) / share_uid / f"{share_uid}.{file_extension}"
    status_file = Path(user_upload_file_base) / share_uid / 'status.json'

    if not user:
        result['message'] = 'Only logged in users can upload datasets.'
        return result

    # formats can be h5ad, rdata, excel, or mex_3tab
    if dataset_format == 'mex_3tab':
        if not filename.endswith('tar.gz') and not filename.endswith('zip'):
            result['message'] = 'Invalid file extension for MEX 3-tab format. Expected .tar.gz or .zip'
            return result

    if dataset_format == 'excel':
        if not filename.endswith('xlsx') and not filename.endswith('xls'):
            result['message'] = 'Invalid file extension for Excel format. Expected .xlsx or .xls'
            return result

    if dataset_format == "h5ad":
        if not filename.endswith('h5ad'):
            result['message'] = 'Invalid file extension for H5AD format. Expected .h5ad'
            return result

    if dataset_format == 'spatial':
        if not filename.endswith('tar.gz'):
            result['message'] = 'Invalid file extension for Spatial format. Expected .tar or .tar.gz'
            return result

        from gear.spatialhandler import SPATIALTYPE2CLASS
        if spatial_format not in SPATIALTYPE2CLASS:
            result['message'] = 'Invalid spatial format specified.'
            return result


    try:
        with open(dataset_filename, 'wb') as f:
            f.write(form['dataset_file'].file.read())
        result['success'] = 1
        result['message'] = 'Dataset file saved successfully.'

        status = {
            "process_id": None,
            "status": "uploaded",
            "message": "The dataset has been uploaded and is pending processing",
            "progress": 0
        }

        with open(status_file, 'w') as f:
            f.write(json.dumps(status))

    except Exception as e:
        result['message'] = 'Error saving dataset file: ' + str(e)

    return result

if __name__ == '__main__':
    result = main()
    print(json.dumps(result))
