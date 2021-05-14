#!/opt/bin/python3
'''
This script checks if a dataset is "Public" or "Private"
If "Public", the data file will be downloaded for the user in .tab format
If "Private", the user will be redirected get the main gEAR page: http://gear.igs.umaryland.edu
'''

from shutil import copyfileobj
import sys
import cgi, json
import os

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    form = cgi.FieldStorage()
    dataset_id = cgi.escape(form.getvalue('dataset_id'))

    file_types = ['.tab', '.csv', '.txt']
    for type in file_types:
        filename = dataset_id + type
        if os.path.isfile('../datasets_uploaded/'+ filename):
            break

    result = { "success": 0 }

    access_level = get_dataset_access_level(cursor, dataset_id)

    if access_level == 'Public':
        result["success"] = 1

        cursor.close()
        cnx.close()

        print("Content-type: application/octet-stream")
        print("Content-Disposition: attachment; filename={0}".format(filename))
        print()
        sys.stdout.flush()

        with open('../datasets_uploaded/' + filename,'rb') as binfile:
            copyfileobj(binfile, sys.stdout.buffer)

    else:
        result["error"] = "Dataset is the private. Unable to download data file."

        cursor.close()
        cnx.close()

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
        """.format(result['error'], 'http://gear.igs.umaryland.edu'))



def get_dataset_access_level(cursor, dataset_id):
    query="""
        SELECT is_public
        FROM dataset
        WHERE id = %s
        """
    cursor.execute(query, (dataset_id,))
    access = None
    for (is_public, ) in cursor:
        if is_public == 0:
            access = 'Private'
        else:
            access = 'Public'
    return access

if __name__ == '__main__':
    main()
