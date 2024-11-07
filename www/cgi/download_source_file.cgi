#!/opt/bin/python3
'''
Given a dataset's ID, this allows for the download of that dataset's tarball
or H5AD file.

'''

import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb


def download_file(file_path, file_name):
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

def main():
    form = cgi.FieldStorage()
    dataset_id = html.escape(form.getvalue('dataset_id'))
    analysis_id = html.escape(form.getvalue('analysis_id', ""))
    session_id = html.escape(form.getvalue('session_id', ""))
    dtype = html.escape(form.getvalue('type'))
    dataset = geardb.Dataset(id=dataset_id)
    tarball_path = dataset.get_tarball_path()
    h5ad_path = dataset.get_file_path()

    # if analysis ID is passed, retrieve the h5ad file for the analysis to download
    if analysis_id:

        # Need session id to get "user_unsaved" analyses
        if not session_id:
          session_id = None

        analysis = geardb.Analysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id)
        analysis.discover_type()
        try:
          h5ad_path = analysis.dataset_path()
        except Exception as e:
          print(str(e), file=sys.stderr)
          h5ad_path = ""

    if dtype == 'tarball' and os.path.isfile(tarball_path):
        download_file(tarball_path, f"{dataset_id}.tar.gz")
    elif dtype == 'h5ad' and os.path.isfile(h5ad_path):
        download_file(h5ad_path, f"{dataset_id}.h5ad")
    else:
        raise FileNotFoundError("File not found")

if __name__ == '__main__':
    main()
