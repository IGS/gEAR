#!/opt/bin/python3

"""
This reads an image file created by a gEAR plotting function and prints it as a binary
stream, suitable for serving as an img['src'] from an html page.
"""

import cgi
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')
    user = geardb.get_user_from_session_id(session_id)

    ana = geardb.Analysis(id=analysis_id, type=analysis_type, dataset_id=dataset_id, 
                          session_id=session_id, user_id=user.id )

    ## if the analysis is primary images are getting saved as user_unsaved
    if ana.type == 'primary':
        ana.type = 'user_unsaved'

    # this has to be what the scanpy library calls it, since we don't seem to be able to
    #  control it: scanpy GitHub ticket #73
    analysis_name = form.getvalue('analysis_name')

    data_file_path = ana.dataset_path()
    ana_directory = os.path.dirname(data_file_path)
    image_path = "{0}/figures/{1}.png".format(ana_directory, analysis_name)
    #print("DEBUG: streaming image at path: {0}".format(image_path), file=sys.stderr)

    with open(image_path, 'rb') as f:
        print("Content-Type: image/png\n")
        sys.stdout.flush() # <---
        sys.stdout.buffer.write(f.read())


if __name__ == '__main__':
    main()
