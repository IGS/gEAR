#!/opt/bin/python3

"""
Used by the expression uploader, this stores the actual dataset from the form
and saves it to a file for processing.

Writes a file at: ../uploads/files/<session_id>/<share_uid>/<share_uid>.<ext>
"""

import cgi
import json
import shutil
import sys
from pathlib import Path

lib_path = Path(__file__).resolve().parents[2] / 'lib'
sys.path.append(str(lib_path))
import geardb
from werkzeug.utils import secure_filename

def main():
    print('Content-Type: application/json\n\n')
    form = cgi.FieldStorage()
    session_id = secure_filename(form.getfirst('session_id', ''))
    share_uid = secure_filename(form.getfirst('share_uid', ''))
    dataset_format = form.getfirst('dataset_format')
    spatial_format = form.getfirst('spatial_format')  # may be None
    # Sent by the browser from File.size, so we can confirm the whole file actually arrived.
    # A dropped/truncated connection can otherwise leave a silently-truncated file on disk with
    # no error surfaced anywhere in the pipeline.
    expected_size = form.getfirst('expected_size')

    if not share_uid: # should never happen
        error_msg = f"Unexpected missing share_uid in store_expression_dataset.cgi. session_id={session_id!r}"
        print(error_msg, file=sys.stderr)
        result = {'success': 0, 'message': 'Internal error: share_uid missing (this should never happen). Please contact support.'}
        return result

    user = geardb.get_user_from_session_id(session_id)
    result = {'success': 0, 'message': ''}

    filename = form['dataset_file'].filename

    # Lowercase the extension so later steps (processing, finalize) can find the file
    #  by a fixed name regardless of how the user's file was named (e.g. .RDS, .XLSX)
    lower_filename = filename.lower()

    if lower_filename.endswith('.tar.gz'):
        file_extension = 'tar.gz'
    else:
        file_extension = secure_filename(lower_filename.split('.')[-1])

    if not file_extension:
        result['message'] = 'Invalid dataset file name.'
        return result

    # This should already have been created when the metadata was stored
    user_upload_file_base = "../uploads/files/{0}".format(session_id)

    dataset_filename = (Path(user_upload_file_base) / share_uid / f"{share_uid}.{file_extension}").resolve()
    status_file = Path(user_upload_file_base) / share_uid / 'status.json'

    uploads_base = Path(user_upload_file_base).resolve()
    if not dataset_filename.is_relative_to(uploads_base):
        result['message'] = 'Invalid dataset file name.'
        return result

    if not user:
        result['message'] = 'Only logged in users can upload datasets.'
        return result

    # formats can be h5ad, rdata, excel, or mex_3tab
    if dataset_format == 'mex_3tab':
        if file_extension not in ('tar.gz', 'zip'):
            result['message'] = 'Invalid file extension for MEX 3-tab format. Expected .tar.gz or .zip'
            return result

    if dataset_format == 'excel':
        if file_extension == 'xls':
            result['message'] = 'Legacy .xls files are not supported. Please re-save the file as .xlsx and upload again.'
            return result
        if file_extension != 'xlsx':
            result['message'] = 'Invalid file extension for Excel format. Expected .xlsx'
            return result

    if dataset_format == "h5ad":
        if file_extension != 'h5ad':
            result['message'] = 'Invalid file extension for H5AD format. Expected .h5ad'
            return result

    if dataset_format == "rds":
        if file_extension != 'rds':
            result['message'] = 'Invalid file extension for RDS format. Expected .rds, .Rds, or .RDS'
            return result

    if dataset_format == 'spatial':
        if file_extension != 'tar.gz':
            result['message'] = 'Invalid file extension for Spatial format. Expected .tar.gz'
            return result

        from gear.spatialhandler import SPATIALTYPE2CLASS
        if spatial_format not in SPATIALTYPE2CLASS:
            result['message'] = 'Invalid spatial format specified.'
            return result


    try:
        # Stream rather than read-then-write -- form['dataset_file'].file is already a real,
        # disk-backed file object at this point (cgi.FieldStorage spools uploads to a temp file),
        # so this avoids needlessly buffering a multi-GB upload in memory a second time.
        with open(dataset_filename, 'wb') as f:
            shutil.copyfileobj(form['dataset_file'].file, f)

        # Confirm the whole file actually arrived. cgi.FieldStorage's multipart parser silently
        # tolerates a dropped/truncated connection (no exception raised), so without this check a
        # partial upload would be reported as a success and processed as if it were complete.
        actual_size = dataset_filename.stat().st_size
        if expected_size is not None and str(actual_size) != str(expected_size):
            dataset_filename.unlink(missing_ok=True)
            result["success"] = 0
            result['message'] = (
                'Upload appears to be incomplete (expected {} bytes, received {}). '
                'This can happen if the connection was interrupted. Please try uploading again.'
            ).format(expected_size, actual_size)
            return result

        result['success'] = 1
        result['message'] = 'Dataset file saved successfully.'

        status = {
            "job_id": None,
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
