#!/opt/bin/python3
'''
Given a dataset's ID, this allows for the download of that dataset's tarball
or H5AD file.

'''

from shutil import copyfileobj
import sys
import cgi, html, json
import os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

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
        print("Content-type: application/octet-stream")
        print("Content-Disposition: attachment; filename={0}.tar.gz".format(dataset_id))
        print()
        sys.stdout.flush()

        with open(tarball_path, 'rb') as binfile:
            copyfileobj(binfile, sys.stdout.buffer)

    elif dtype == 'h5ad' and os.path.isfile(h5ad_path):
        print("Content-type: application/octet-stream")
        print("Content-Disposition: attachment; filename={0}.h5ad".format(dataset_id))
        print()
        sys.stdout.flush()

        with open(h5ad_path, 'rb') as binfile:
            copyfileobj(binfile, sys.stdout.buffer)

    else:
        raise FileNotFoundError("File not found")

if __name__ == '__main__':
    main()
