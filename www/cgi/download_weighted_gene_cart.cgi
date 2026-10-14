#!/opt/bin/python3

"""
Downloads a weighted gene cart
"""

import cgi
import sys
from shutil import copyfileobj
from pathlib import Path

from werkzeug.utils import secure_filename

lib_path = Path(__file__).resolve().parents[2].joinpath('lib')
sys.path.insert(0, str(lib_path))

import geardb

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

def main():

    form = cgi.FieldStorage()
    share_id = form.getfirst('share_id')
    share_id = secure_filename(share_id or "")

    if not share_id:
        raise ValueError("Share ID not provided")

    # Get the gene symbols from the shared cart file
    file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format("cart." + share_id))
    if not str(file_path).startswith(str(CARTS_BASE_DIR)):
        raise ValueError("Not allowed.")

    # Check before printing headers, so a missing file can still get a 404
    if not file_path.is_file():
        raise FileNotFoundError(f"No weighted gene list file found for share ID {share_id}")

    print("Content-type: application/octet-stream")
    # Only the file name; this used to send the full server-side path
    print(f"Content-Disposition: attachment; filename={file_path.name}")
    print()
    sys.stdout.flush()

    with open(file_path, 'rb') as binfile:
        copyfileobj(binfile, sys.stdout.buffer)

def print_error(status, message):
    """Print a plain-text error response with the given HTTP status (before any other headers)."""
    print(f"Status: {status}")
    print("Content-Type: text/plain")
    print()
    print(message)

if __name__ == '__main__':
    try:
        main()
    except FileNotFoundError as e:
        print_error("404 Not Found", str(e))
    except ValueError as e:
        print_error("400 Bad Request", str(e))
