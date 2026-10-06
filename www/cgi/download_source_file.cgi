#!/opt/bin/python3
'''
Given a dataset's ID, this allows for the download of that dataset's tarball
or H5AD file.

The tarball and H5AD (including analysis H5ADs) are only served when the dataset
is marked downloadable or the requester owns it. The requester is identified from
the session_id parameter or, for plain download links, the gear_session_id cookie.
Errors return a plain-text message with HTTP 400, 403 or 404.
'''

import cgi
import os
import sys
from http.cookies import SimpleCookie

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.analysis import Analysis
from werkzeug.utils import secure_filename

def download_file(file_path, file_name):
    """
    Stream a file to stdout as an attachment download in 8KB chunks.
    """
    print("Content-type: application/octet-stream")
    print(f"Content-Disposition: attachment; filename={file_name}")
    print()
    sys.stdout.flush()

    with open(file_path, 'rb') as binfile:
        while True:
            chunk = binfile.read(8192)  # Read in 8KB chunks
            if not chunk:
                break
            sys.stdout.buffer.write(chunk)
            sys.stdout.buffer.flush()

def to_file(content, prefix='', suffix=''):
    """
    Write string content to a named temporary file.

    Returns:
        Path to the temporary file (caller is responsible for deleting it).
    """
    import tempfile
    temp = tempfile.NamedTemporaryFile(delete=False, prefix=prefix, suffix=suffix, mode='w')
    temp.write(content)
    temp.close()
    return temp.name

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

def print_error(status, message):
    """Print a plain-text error response with the given HTTP status (before any other headers)."""
    print(f"Status: {status}")
    print("Content-Type: text/plain")
    print()
    print(message)

def main():
    form = cgi.FieldStorage()
    dataset_id = secure_filename(form.getfirst('dataset_id') or "")
    share_id = secure_filename(form.getfirst("share_id") or "")
    analysis_id = secure_filename(form.getfirst('analysis_id') or "")
    session_id = get_request_session_id(form)
    dtype = form.getfirst('type') or ""

    if not dataset_id and not share_id:
        raise ValueError("Either dataset ID or share ID must be provided")

    if not dtype:
        raise ValueError("Type must be provided (tarball, h5ad, or metadata)")
    elif dtype == "metadata":
        metadata_content = geardb.get_metadata_by_share_id(share_id)
        if metadata_content:
            temp_file_path = to_file(metadata_content, prefix=f"{share_id}_metadata", suffix=".csv")
            try:
                download_file(temp_file_path, f"{share_id}.metadata.csv")
            finally:
                os.remove(temp_file_path)
        else:
            raise FileNotFoundError("Metadata not found")
    else:
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
            user = geardb.get_user_from_session_id(session_id) if session_id else None
            if user is None or user.id != dataset.owner_id:
                raise PermissionError("This dataset is not available for download.")

        archive_path = dataset.get_source_archive_path()
        h5ad_path = dataset.get_file_path()

        # if analysis ID is passed, retrieve the h5ad file for the analysis to download
        if analysis_id:
            # Need session id to get "user_unsaved" analyses
            if not session_id:
                session_id = None

            analysis = Analysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id)
            analysis.discover_type()
            try:
                h5ad_path = analysis.dataset_path
            except Exception as e:
                print(str(e), file=sys.stderr)
                h5ad_path = ""

        if dtype == 'tarball' and archive_path:
            # Keep the archive's own extension (.tar.gz, .tar or .zip) in the download name
            archive_extension = ".tar.gz" if archive_path.endswith(".tar.gz") else os.path.splitext(archive_path)[1]
            download_file(archive_path, f"{share_id}{archive_extension}")
        elif dtype == 'h5ad' and os.path.isfile(h5ad_path):
            download_file(h5ad_path, f"{share_id}.h5ad")
        else:
            raise FileNotFoundError("File not found")

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