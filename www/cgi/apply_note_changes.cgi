#!/opt/bin/python3

"""
This script handles notes by:
1) Saving changes to existing notes
2) Creating new notes
3) Removing notes from the database

The action taken is determined by the 'scope' parameter, which can be:
'new' for creating new notes, or
'edit' for saving changes to existing notes
'remove' for deleting notes from the database

"""

import cgi, json
from datetime import datetime
import os,sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')
    result = {'success': 0, 'notes': [] }

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    title = form.getvalue('title')
    ldesc = form.getvalue('ldesc')
    dataset_id = form.getvalue('dataset_id') # 'None' if scope='edit'
    is_public = form.getvalue('access_level')
    scope = form.getvalue('scope') # 'new', 'edit', or 'remove'
    note_id = form.getvalue('note_id') #if scope = 'edit', this is contains the note_id

    current_user_id = get_user_id_from_session_id(cursor, session_id)

    # Is a user even logged in?
    if current_user_id is None:
        result['error'] = 'Please log in to access or create notes for this dataset.'

    else:
        # A user is logged in

        # Add the note to the database
        if scope == 'new':
            add_new_note(cursor, title, ldesc, current_user_id, dataset_id, is_public)
            cnx.commit()
            result['success'] = 1

        # Save changes to existing note
        if scope == 'edit':
            if note_id is not None:
                save_changes(cursor, note_id, title, ldesc, is_public)
                cnx.commit()
                result['success'] = 1
            else:
                result['error'] = 'Not able to save changes to note. Invalid note ID.'

        # Remove the note from database
        if scope == 'remove':
            #Does user own the note?
            is_owner = check_note_ownership(cursor, note_id, current_user_id)
            if is_owner == True:
                remove_note(cursor, note_id)
                cnx.commit()
                result['success'] = 1
            else:
                result['error'] = 'Not able to remove note. Not note owner.'


    cursor.close()
    cnx.close()

    print(json.dumps(result))

def remove_note(cursor, note_id):
    qry = """
        DELETE FROM note
        WHERE id = %s
    """
    cursor.execute(qry, (note_id,))

def check_note_ownership(cursor, note_id, current_user_id):
    qry = """
        SELECT user_id
        FROM note
        WHERE id = %s
    """
    cursor.execute(qry, (note_id,))
    for (user_id,) in cursor:
        if user_id == current_user_id:
            return True
        else:
            return False

def save_changes(cursor, note_id, title, ldesc, is_public):
    qry = """
       UPDATE note
       SET title = %s, ldesc = %s, is_public = %s, date_last_changed = NOW()
       WHERE id = %s
    """
    cursor.execute(qry, (title, ldesc, is_public, note_id,))

def add_new_note(cursor, title, ldesc, current_user_id, dataset_id, is_public):
    qry = """
       INSERT INTO note (title, ldesc, user_id, dataset_id, is_public, date_added, date_last_changed)
       VALUES (%s, %s, %s, %s, %s, NOW(), NOW() )
    """
    cursor.execute(qry, (title, ldesc, current_user_id, dataset_id, is_public,))


def get_user_id_from_session_id(cursor, session_id):
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid

    return user_id

if __name__ == '__main__':
    main()
