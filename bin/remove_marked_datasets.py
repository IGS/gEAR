#!/opt/bin/python3

"""
Removes any dataset from gEAR database where dataset.marked_for_removal = 1

"""

import cgi, json
import configparser
import mysql.connector
from mysql.connector import errorcode
import sys
import os

def main():
    config = configparser.ConfigParser()
    config.read('../gear.ini')

    cnx = None
    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'])
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password", file=sys.stderr)
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist", file=sys.stderr)
        else:
            print(err, file=sys.stderr)

    print("\nConnected to gEAR database.\n")

    cursor = cnx.cursor()
    marked_dataset_ids = list()
    query = """
        SELECT id
        FROM dataset
        WHERE marked_for_removal = 1
    """

    cursor.execute(query,)
    for row in cursor:
        marked_dataset_ids.append(row[0])

    dataset_count = len(marked_dataset_ids)


    if dataset_count == 0:
        print("No datasets identified for removal.\n")
        sys.exit()

    else:
        print( "{0} datasets identified for removal.\n".format(dataset_count) )

        # Make list of files in datasets_uploaded
        datasets_uploaded = os.listdir('./www/datasets_uploaded')

        # Delete each dataset
        for dataset in marked_dataset_ids:
            print( "Removing dataset: {0}...".format(dataset) )
            remove_from_expression(cnx, cursor, dataset)
            remove_from_layout_members(cursor, dataset)
            remove_from_dataset_shares(cursor, dataset)
            remove_from_dataset_tag(cursor, dataset)
            remove_from_note(cursor, dataset)
            remove_from_dataset(cursor, dataset)

            cnx.commit()

            file_uploaded_count = 0
            for file in datasets_uploaded:
                # Remove the dataset's uploaded files
                if dataset in file:
                    os.remove('./www/datasets_uploaded/' + file)
                    file_uploaded_count += 1
            print("...{0} files removed from 'datasets_uploaded' directory.".format(file_uploaded_count))

            dataset_count -= 1
            print( "Number of datasets left to remove: {0}\n".format(dataset_count) )

        cursor.close()
        cnx.close()

        print("Dataset removal complete.")


def remove_from_expression(conn, cursor, dataset_id):
    qry = """
        DELETE FROM expression
        WHERE dataset_id = %s
        LIMIT 10000
    """
    while True:
        cursor.execute(qry, (dataset_id,))
        conn.commit()
        print("... deleted 10000 rows from expression for dataset {0}".format(dataset_id))
        if cursor.rowcount < 10000:
            break

    print("...'expression' data removed.")

def remove_from_layout_members(cursor, dataset_id):
    qry = """
        DELETE FROM layout_members
        WHERE dataset_id = %s
    """
    cursor.execute(qry, (dataset_id,))
    print("...'layout_members' data removed.")

def remove_from_dataset_shares(cursor, dataset_id):
    qry = """
        DELETE FROM dataset_shares
        WHERE dataset_id = %s
    """
    cursor.execute(qry, (dataset_id,))
    print("...'dataset_shares' data removed.")

def remove_from_dataset_tag(cursor, dataset_id):
    qry = """
        DELETE FROM dataset_tag
        WHERE dataset_id = %s
    """
    cursor.execute(qry, (dataset_id,))
    print("...'dataset_tag' data removed.")

def remove_from_note(cursor, dataset_id):
    qry = """
        DELETE FROM note
        WHERE dataset_id = %s
    """
    cursor.execute(qry, (dataset_id,))
    print("...'note' data removed.")

def remove_from_dataset(cursor, dataset_id):
    qry = """
        DELETE FROM dataset
        WHERE id = %s
    """
    cursor.execute(qry, (dataset_id,))
    print("...'dataset' data removed.")


if __name__ == '__main__':
    main()
