#!/opt/bin/python3

"""
Deletes a comment.
"""

import cgi
import configparser
import json
import math
import mysql.connector
from mysql.connector import errorcode
import re
import sys
from xml.dom import minidom
import datetime

sys.path.append('../..')
import loaderutils

def main():
    print('Content-Type: application/json\n\n')

    config = configparser.ConfigParser()
    config.read('../../gear.ini')
    result = { 'success':0 }

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

    cursor = cnx.cursor()

    form = cgi.FieldStorage()
    comment_id = form.getvalue('comment_id')

    result = { 'success': 0 }

    #delete comment from comment_tag - Do this first to avoid constraint fail
    delete_from_comment_tag(cursor, comment_id)

    #delete comment
    delete_comment(cursor, comment_id)
    cnx.commit()

    cursor.close()
    cnx.close()

    result['success'] = 1

    print(json.dumps(result))


def delete_from_comment_tag(cursor, comment_id):
    qry = """
        DELETE FROM comment_tag
        WHERE comment_id = %s
    """
    cursor.execute(qry, (comment_id,))


def delete_comment(cursor, comment_id):
    qry = """
        DELETE FROM comment
        WHERE id = %s
    """
    cursor.execute(qry, (comment_id,))



if __name__ == '__main__':
    main()
