#!/opt/bin/python3

"""
Generate help_ids for guser entries where help_id == NULL.

"""

import cgi, json
import configparser
import mysql.connector
from mysql.connector import errorcode
import sys
import uuid

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
        FROM guser
        WHERE help_id is NULL
    """

    cursor.execute(query,)

    guser_ids=list()
    for row in cursor:
        guser_ids.append(row[0])

    print("{0} help_ids needed.".format(len(guser_ids)))

    for user_id in guser_ids:
        help_id = str(uuid.uuid4())
        print(help_id)
        qry = """
            UPDATE guser
            SET help_id = %s
            WHERE id = %s
        """
        cursor.execute(qry, (help_id, user_id,))
        cnx.commit()
        print("help_id added to guser: {0}".format(user_id))


    cursor.close()
    cnx.close()
    print("Task Complete. help_ids generated and inserted into guser table.")


if __name__ == '__main__':
    main()
