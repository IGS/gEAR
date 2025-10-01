#!/opt/bin/python3
'''
Given a dataset's ID, this allows for the download of that dataset's tarball
or H5AD file.

'''

import cgi
import html
import io
import sys
import zipfile
from pathlib import Path

lib_path = Path(__file__).resolve().parents[2].joinpath('lib')
sys.path.append(str(lib_path))
import geardb

PROJECTION_DATASET_DIR = (Path(__file__).resolve().parent.parent / 'projections' / 'by_dataset').resolve()

def main():
    form = cgi.FieldStorage()
    dataset_id = html.escape(form.getvalue('dataset_id', ""))
    share_id = html.escape(form.getvalue("share_id", ""))
    projection_id = html.escape(form.getvalue('projection_id', ""))

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

    coeff_path = PROJECTION_DATASET_DIR / dataset_id / f"{projection_id}.csv"
    pval_path = PROJECTION_DATASET_DIR / dataset_id / f"{projection_id}_pval.csv"

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



if __name__ == '__main__':
    main()
