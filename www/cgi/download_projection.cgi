#!/opt/bin/python3
'''
Given a dataset and projection ID, downloads the projection output (coefficients,
plus p-values when present) as a zip file.

Like download_source_file.cgi, the output is only served when the dataset is marked
downloadable or the requester owns it (session_id parameter or gear_session_id cookie).
Errors return a plain-text message with HTTP 400, 403 or 404.
'''

import cgi
import io
import os
import sys
import zipfile
from http.cookies import SimpleCookie
from pathlib import Path

lib_path = Path(__file__).resolve().parents[2].joinpath('lib')
sys.path.append(str(lib_path))
import geardb
from werkzeug.utils import secure_filename

PROJECTION_DATASET_DIR = (Path(__file__).resolve().parent.parent / 'projections' / 'by_dataset').resolve()

def get_request_session_id(form):
    """
    Return the requester's session ID from the session_id parameter, falling back to the
    gear_session_id cookie (download links opened in the browser don't pass session_id).
    """
    session_id = form.getfirst('session_id') or ""
    if not session_id:
        cookie = SimpleCookie(os.environ.get("HTTP_COOKIE", ""))
        if "gear_session_id" in cookie:
            session_id = cookie["gear_session_id"].value
    return secure_filename(session_id)

def main():
    form = cgi.FieldStorage()
    dataset_id = secure_filename(form.getfirst('dataset_id', ""))
    share_id = secure_filename(form.getfirst("share_id", ""))
    projection_id = secure_filename(form.getfirst('projection_id', ""))

    if not dataset_id and not share_id:
        raise ValueError("Either dataset ID or share ID must be provided")

    # if share ID is passed, retrieve the dataset by share ID
    if share_id:
        dataset = geardb.get_dataset_by_share_id(share_id, False)
        if not dataset:
            raise FileNotFoundError(f"Dataset not found for the provided share ID {share_id}")
        dataset_id = dataset.id

    elif dataset_id:
        # Legacy support
        dataset = geardb.get_dataset_by_id(dataset_id, False)
        if not dataset:
            raise FileNotFoundError(f"Dataset not found for the provided dataset ID {dataset_id}")
        share_id = dataset.share_id

    # Honor the dataset's "Is downloadable" setting; the owner can always download
    if not dataset.is_downloadable:
        session_id = get_request_session_id(form)
        user = geardb.get_user_from_session_id(session_id) if session_id else None
        if user is None or user.id != dataset.owner_id:
            raise PermissionError("This dataset is not available for download.")

    coeff_path = (PROJECTION_DATASET_DIR / dataset_id / f"{projection_id}.csv").resolve()
    pval_path = (PROJECTION_DATASET_DIR / dataset_id / f"{projection_id}_pval.csv").resolve()

    if not coeff_path.is_relative_to(PROJECTION_DATASET_DIR) or not pval_path.is_relative_to(PROJECTION_DATASET_DIR):
        raise ValueError("Invalid dataset ID or projection ID.")

    zip_buffer = io.BytesIO()

    # Place coeff file and pval file (if found) in a zip and send that
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        if Path(coeff_path).is_file():
            zf.write(coeff_path, f"{projection_id}.coeff.csv")
        else:
            raise FileNotFoundError("Projection output not found")

        # This file is optional. If added, zip download becomes a directory instead of a single file
        if Path(pval_path).is_file():
            zf.write(pval_path, f"{projection_id}.pval.csv")

    # Rewind the buffer
    zip_buffer.seek(0)

    # Download the zip
    print("Content-type: application/octet-stream")
    print(f"Content-Disposition: attachment; filename={projection_id}.zip")
    print()
    sys.stdout.flush()
    # Stream the buffer to stdout
    while True:
        chunk = zip_buffer.read(8192)
        if not chunk:
            break
        sys.stdout.buffer.write(chunk)
        sys.stdout.buffer.flush()



def print_error(status, message):
    """Print a plain-text error response with the given HTTP status (before any other headers)."""
    print(f"Status: {status}")
    print("Content-Type: text/plain")
    print()
    print(message)

if __name__ == '__main__':
    # These are raised before any headers are printed, so a proper status can still be sent
    try:
        main()
    except PermissionError as e:
        print_error("403 Forbidden", str(e))
    except FileNotFoundError as e:
        print_error("404 Not Found", str(e))
    except ValueError as e:
        print_error("400 Bad Request", str(e))
