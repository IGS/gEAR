#!/opt/bin/python3
'''
Given a dataset's ID, this allows for the download of that dataset's tarball
or H5AD file.

'''

from shutil import copyfileobj
import sys
import cgi, json
import os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    dataset_id = cgi.escape(form.getvalue('dataset_id'))
    dtype = cgi.escape(form.getvalue('type'))
    dataset = geardb.Dataset(id=dataset_id)
    tarball_path = dataset.get_tarball_path()
    h5ad_path = dataset.get_file_path()

    if dtype == 'tarball' and os.path.isfile(tarball_path):
        print("Content-type: application/octet-stream")
        print("Content-Disposition: attachment; filename={0}.tar.gz".format(dataset_id))
        print()
        sys.stdout.flush()

        with open(tarball_path, 'rb') as binfile:
            copyfileobj(binfile, sys.stdout.buffer)

    elif dtype == 'h5ad' and os.path.isfile(tarball_path):
        print("Content-type: application/octet-stream")
        print("Content-Disposition: attachment; filename={0}.h5ad".format(dataset_id))
        print()
        sys.stdout.flush()

        with open(h5ad_path, 'rb') as binfile:
            copyfileobj(binfile, sys.stdout.buffer)

    else:
        result_error = "Dataset tarball could not be found. Unable to download data file."

        print("Content-type: text/html")
        print()
        print("""
        <!DOCTYPE html>
        <html>
          <head>
            <meta http-equiv='refresh' content='5;url=http://gear.igs.umaryland.edu/'>
          </head>
          <body>
            <p>Error: {0}</p>
            <p>Redirecting... <a href='{1}'>Click here if you are not redirected</a>
          </body>
        </html>
        """.format(result_error, 'http://gear.igs.umaryland.edu'))

if __name__ == '__main__':
    main()
