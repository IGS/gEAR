#!/opt/bin/python3

"""
Marks a comment as read (1).


mysql> SELECT comment.*, GROUP_CONCAT(tag.label)
FROM comment
JOIN comment_tag ON comment.id = comment_tag.comment_id
JOIN tag ON tag.id = comment_tag.tag_id
GROUP BY comment.id;
+----+------------+-----------+---------+-------------------+------------+---------------------------+---------+---------------------+-------------------------+
| id | first_name | last_name | user_id | email             | title      | message                   | is_read | date_added          | GROUP_CONCAT(tag.label) |
+----+------------+-----------+---------+-------------------+------------+---------------------------+---------+---------------------+-------------------------+
|  3 | asfasdf    | asdfasdfa |    NULL | asdf@gmail.com    | Test       | Test first                |       0 | 2016-06-28 11:48:30 | mouse                   |
|  9 | asfasdf    | asdfasdfa |    NULL | asdf@gmail.com    | Test       | Test first                |       0 | 2016-06-28 12:19:46 | mouse                   |
| 10 | asfasdf    | asdfasdfa |    NULL | asdf@gmail.com    | Test       | Test first                |       0 | 2016-06-28 12:19:48 | mouse                   |
| 11 | asfasdf    | asdfasdfa |    NULL | asdf@gmail.com    | Test       | Test first                |       0 | 2016-06-28 12:19:58 | mouse                   |
| 12 | asdf       | asdf      |    NULL | asdf@gmail.com    | Test       | testing                   |       0 | 2016-06-28 12:20:30 | mouse                   |
| 13 | asdf       | asdf      |    NULL | asdf@gmail.com    | Test       | testing                   |       0 | 2016-06-28 12:20:49 | mouse                   |
| 15 | Dustin     | Olley     |    NULL | d@gmail.com       | Test       | test                      |       0 | 2016-06-28 12:22:28 | sox2                    |
| 16 | asd        | asd       |    NULL | ass               | as         | as                        |       0 | 2016-06-28 12:47:46 | mouse                   |
| 17 | asdf       | asfd      |    NULL | asdf              | asdf       | asdf                      |       0 | 2016-06-28 13:20:10 | mouse                   |
| 18 | Dustin     | Olley     |    NULL | hi@hello.com      | Multi Tags | tags tags                 |       0 | 2016-06-28 13:21:57 | tag1, tag2              |
| 19 | asdf       | asdf      |    NULL | asdf              | asdf       | asdf                      |       0 | 2016-06-28 13:25:40 | tag3,tag4               |
| 20 | Token      | Tag       |    NULL | t.tag@gmail.com   | Token Tag  | Token Tag                 |       0 | 2016-06-29 10:52:43 | spaced tag,newtag,mouse |
| 21 | Token      | Tag       |    NULL | t.tag@gmail.com   | Token Tag  | Token Tag                 |       0 | 2016-06-29 10:52:43 | newtag,mouse,spaced tag |
| 22 | Dustin     | Olley     |    NULL | dolleyj@gmail.com | TEST final | This will work. Then done |       0 | 2016-06-29 11:10:47 | bananas,help,mouse      |
| 23 | Dustin     | Olley     |    NULL | dolleyj@gmail.com | TEST       | This is only a test       |       0 | 2016-08-01 13:30:24 | tag1                    |
+----+------------+-----------+---------+-------------------+------------+---------------------------+---------+---------------------+-------------------------+
Source for putting tags together: http://stackoverflow.com/a/3664419/2900840
"""

import cgi
import json
import math
import re, sys, os

from xml.dom import minidom
import datetime

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
import loaderutils

def main():
    print('Content-Type: application/json\n\n')
    result = { 'success':0, 'comment': list() }

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    form = cgi.FieldStorage()
    comment_id = form.getvalue('comment_id')
    change_is_read = form.getvalue('change_is_read')

    #change comment's read status
    update_comment(cursor, change_is_read, comment_id)
    cnx.commit()

    #get the comment
    result['comment'].extend(get_comment(cursor, comment_id))
    cursor.close()
    cnx.close()
    result['success'] = 1

    print(json.dumps(result))



def update_comment(cursor, change_is_read, comment_id):
    qry = """
        UPDATE comment
        SET comment.is_read = %s
        WHERE comment.id = %s
    """
    cursor.execute(qry, (change_is_read, comment_id))


def get_comment(cursor, comment_id):
    qry = """
        SELECT c.id, c.first_name, c.last_name, c.user_id, c.email, c.title, c.message, c.is_read, c.date_added, GROUP_CONCAT(tag.label)
        FROM comment c
        JOIN comment_tag ON comment_tag.comment_id = c.id
        JOIN tag ON tag.id = comment_tag.tag_id
        WHERE c.id = %s
        """
    cursor.execute(qry, (comment_id,))
    comment = list()

    for row in cursor:
        if row[7] == 0:
            is_read = False
        else:
            is_read = True

        # No tags for entry
        if row[9] is None:
            comment.append({
            'id': row[0],
            'first_name': row[1],
            'last_name': row[2],
            'user_id': row[3],
            'email': row[4],
            'title': row[5],
            'message': row[6],
            'is_read': is_read,
            'date_added': str(row[8]),
            'tags': None
            })
        else:
            comment.append({
            'id': row[0],
            'first_name': row[1],
            'last_name': row[2],
            'user_id': row[3],
            'email': row[4],
            'title': row[5],
            'message': row[6],
            'is_read': is_read,
            'date_added': str(row[8]),
            'tags': row[9].split(',')
            })

    return comment



if __name__ == '__main__':
    main()
