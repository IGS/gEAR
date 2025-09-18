#!/opt/bin/python3

"""
Used by the expression uploader, this stores the actual dataset from the form
and saves it to a file for processing.

Writes a file at: ../uploads/files/<session_id>/<share_uid>/<share_uid>.<ext>
"""

import cgi
import json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_uid = form.getvalue('dataset_uid')
    share_uid = form.getvalue('share_uid')
    dataset_format = form.getvalue('dataset_format')
    spatial_format = form.getvalue('spatial_format')  # may be None

    user = geardb.get_user_from_session_id(session_id)
    result = {'success': 0, 'message': ''}

    filename = form['dataset_file'].filename

    if filename.endswith('.tar.gz'):
        file_extension = 'tar.gz'
    else:
        file_extension = filename.split('.')[-1]

    # This should already have been created when the metadata was stored
    user_upload_file_base = "../uploads/files/{0}".format(session_id)
    dataset_filename = os.path.join(user_upload_file_base, share_uid, share_uid + '.' + file_extension)
    status_file = os.path.join(user_upload_file_base, share_uid, 'status.json')

    if not user:
        result['message'] = 'Only logged in users can upload datasets.'
        print(json.dumps(result))
        sys.exit(0)

    # formats can be h5ad, rdata, excel, or mex_3tab
    if dataset_format == 'mex_3tab':
        if not filename.endswith('tar.gz') and not filename.endswith('zip'):
            result['message'] = 'Invalid file extension for MEX 3-tab format. Expected .tar.gz or .zip'
            print(json.dumps(result))
            sys.exit(0)

    if dataset_format == 'excel':
        if not filename.endswith('xlsx') and not filename.endswith('xls'):
            result['message'] = 'Invalid file extension for Excel format. Expected .xlsx or .xls'
            print(json.dumps(result))
            sys.exit(0)

    if dataset_format == 'spatial':
        result["message"] = "NOT YET IMPLEMENTED"
        print(json.dumps(result))
        sys.exit(0)

        if not filename.endswith('tar.gz'):
            result['message'] = 'Invalid file extension for Spatial format. Expected .tar or .tar.gz'
            print(json.dumps(result))
            sys.exit(0)
        from gear.spatialhandler import SPATIALTYPE2CLASS
        if spatial_format not in SPATIALTYPE2CLASS:
            result['message'] = 'Invalid spatial format specified.'
            print(json.dumps(result))
            sys.exit(0)

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

    print(json.dumps(result))

if __name__ == '__main__':
    main()
