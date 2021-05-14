#!/opt/bin/python3

"""
Deletes a comment.
"""

import cgi
import json
import math
import re, os, sys
from xml.dom import minidom
import datetime

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
import loaderutils

def main():
    print('Content-Type: application/json\n\n')
    result = { 'success':0 }

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
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
