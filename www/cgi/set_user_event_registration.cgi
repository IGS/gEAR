#!/opt/bin/python3

"""
This script can be used to both add and remove a user's registration status
to an event.

Input: session_id, event_id, registration_status (1 = register, 0 = unregister).

Registrations are counted in sign-up order (see get_event_registration_list.cgi): the
first max_attendees registrations have a seat and the next waitlist_size are waitlisted.
Registering when the seats are taken adds the user to the waitlist; registering when the
waitlist is also full is refused.  Registering twice doesn't add a second registration.
"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    event_id = form.getfirst('event_id')
    registration_status = form.getfirst('registration_status')
    user = geardb.get_user_from_session_id(session_id) if session_id else None

    print('Content-Type: application/json\n\n')

    if not user:
        print(json.dumps({'success': 0, 'msg': 'You must be logged in to register or unregister'}))
        sys.exit(0)

    registration_status = (registration_status or '').strip()
    if registration_status not in ('0', '1'):
        print(json.dumps({'success': 0, 'msg': 'Unrecognized registration status'}))
        sys.exit(0)

    try:
        event_id = int(event_id)
    except (TypeError, ValueError):
        event_id = None

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    # Locks the event row until commit/rollback, so concurrent sign-ups for the same
    #  event are counted one at a time and can't oversubscribe it
    event = get_event_info(cursor, event_id) if event_id is not None else None

    if event is None:
        cnx.rollback()
        print(json.dumps({'success': 0, 'msg': 'Unknown event'}))

    elif registration_status == '1':
        current_count = get_current_count(cursor, event_id)

        if is_registered(cursor, event_id, user.id):
            cnx.rollback()
            print(json.dumps({'success': 1, 'current_count': current_count,
                              'event': event, 'msg': 'You are already registered'}))

        elif current_count >= event['max_attendees'] + event['waitlist_size']:
            cnx.rollback()
            print(json.dumps({'success': 0, 'msg': 'This event is full'}))

        else:
            qry = "INSERT INTO event_registration (event_id, user_id) VALUES (%s, %s)"
            cursor.execute(qry, (event_id, user.id))
            cnx.commit()

            if current_count >= event['max_attendees']:
                msg = 'All seats are taken, so you have been added to the waitlist'
            else:
                msg = 'You have successfully registered'

            print(json.dumps({'success': 1, 'current_count': current_count + 1,
                              'event': event, 'msg': msg}))

    else:
        qry = "DELETE FROM event_registration WHERE event_id = %s AND user_id = %s"
        cursor.execute(qry, (event_id, user.id))
        cnx.commit()
        print(json.dumps({'success': 1, 'current_count': get_current_count(cursor, event_id),
                          'event': event, 'msg': 'You have successfully unregistered'}))

    cursor.close()
    cnx.close()

def get_current_count(cursor, event_id):
    """
    Return the number of registrations for an event.
    """
    qry = """
          SELECT COUNT(er.id)
            FROM event_registration er
           WHERE er.event_id = %s
    """
    cursor.execute(qry, (event_id,))
    rows = cursor.fetchall()   # read the whole result (the cursor is unbuffered)

    return rows[0][0] if rows else 0

def is_registered(cursor, event_id, user_id):
    """
    Return True if the user already has a registration (seat or waitlist) for the event.
    """
    qry = "SELECT id FROM event_registration WHERE event_id = %s AND user_id = %s LIMIT 1"
    cursor.execute(qry, (event_id, user_id))
    return len(cursor.fetchall()) > 0

def get_event_info(cursor, event_id):
    """
    Return a dict of event capacity info (max_attendees, waitlist_size) with zeroed counters,
    or None if there is no such event.  Locks the event row for the current transaction.
    """
    qry = "SELECT label, max_attendees, waitlist_size FROM event WHERE id = %s FOR UPDATE"
    cursor.execute(qry, (event_id,))
    rows = cursor.fetchall()

    if not rows:
        return None

    return {'attendees': 0,
            'user_attending': 0,
            'user_waitlisted': 0,
            'max_attendees': rows[0][1],
            'waitlist_size': rows[0][2]}


if __name__ == '__main__':
    main()
