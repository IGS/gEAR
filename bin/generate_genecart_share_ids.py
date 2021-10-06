#!/opt/bin/python3

"""
Generate share IDs for gene carts where they are currently NULL

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
        FROM gene_cart
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
            UPDATE gene_cart
            SET share_id = %s
            WHERE id = %s
        """
        cursor.execute(qry, (share_id, id,))
        cnx.commit()



if __name__ == '__main__':
    main()
