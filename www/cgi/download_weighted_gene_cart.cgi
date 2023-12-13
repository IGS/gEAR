#!/opt/bin/python3

"""
Downloads a weighted gene cart
"""

import cgi
import sys
from shutil import copyfileobj
from pathlib import Path

lib_path = Path(__file__).resolve().parents[2].joinpath('lib')
sys.path.insert(0, str(lib_path))

import geardb

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

def main():

    form = cgi.FieldStorage()
    share_id = form.getvalue('share_id')

    if not share_id:
        raise Exception("ERROR: Share ID not provided")

    # Get the gene symbols from the shared cart file
    file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format("cart." + share_id)).resolve()
    if not str(file_path).startswith(str(CARTS_BASE_DIR)):
        raise ValueError("Not allowed.")

    print("Content-type: application/octet-stream")
    print(f"Content-Disposition: attachment; filename={file_path}")
    print()
    sys.stdout.flush()

    with open(file_path, 'rb') as binfile:
        copyfileobj(binfile, sys.stdout.buffer)

if __name__ == '__main__':
    main()
