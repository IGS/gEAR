#!/opt/bin/python3

"""
Generate share IDs for profiles where they are currently NULL

"""

import cgi
import json
import os
import sys
import uuid

lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    cnx = geardb.Connection()

    cursor = cnx.get_cursor()
    marked_dataset_ids = list()
    query = """
        SELECT id
        FROM layout
        WHERE share_id is NULL
    """

    cursor.execute(query,)

    ids=list()
    for row in cursor:
        ids.append(row[0])

    for id in ids:
        # Get just the first block of the UUID
        share_id = str(uuid.uuid4()).split('-')[0]
        print(share_id)
        qry = """
            UPDATE layout
            SET share_id = %s
            WHERE id = %s
        """
        cursor.execute(qry, (share_id, id,))
        cnx.commit()


    #cursor.close()
    #cnx.close()


if __name__ == '__main__':
    main()
