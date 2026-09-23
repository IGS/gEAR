#!/opt/bin/python3

"""
get_event_registration_list.cgi - Report event registration status for the user.

Input: session_id (required), min_event_id, max_event_id.
Output: JSON {event_id: {attendees, max_attendees, waitlist_size, user_attending,
        user_waitlisted}}, or {} if no valid session.
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    min_event_id = form.getvalue('min_event_id')
    max_event_id = form.getvalue('max_event_id')
    user = geardb.get_user_from_session_id(session_id) if session_id else None

    print('Content-Type: application/json\n\n')

    if not user:
        print("{}")
        sys.exit(0)

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    events = get_event_info(cursor, min_event_id, max_event_id)

    qry = """
          SELECT e.id, er.id, u.id
            FROM event e
                 JOIN event_registration er ON er.event_id=e.id
                 JOIN guser u on er.user_id=u.id
           WHERE e.id BETWEEN %s AND %s
    """

    qry_args = [min_event_id, max_event_id]
    
    cursor.execute(qry, qry_args)

    for row in cursor:
        registered_user_id = row[2]
        event = events[row[0]]
        event['attendees'] += 1

        ## is it the current user?
        if registered_user_id == user.id:
            # else are we already over max & waitlist
            if event['attendees'] >= event['max_attendees']:
                event['user_waitlisted'] = 1
            else:
                event['user_attending'] = 1

    print(json.dumps(events))

    cursor.close()
    cnx.close()
    
def get_event_info(cursor, min_event_id, max_event_id):
    """
    Return base info for events with IDs in the given range.

    Returns:
        Dict keyed by event ID with attendee counters initialized to 0.
    """
    qry = "SELECT id, label, max_attendees, waitlist_size FROM event WHERE id BETWEEN %s AND %s"
    events = dict()

    cursor.execute(qry, (min_event_id, max_event_id))

    for row in cursor:
        events[row[0]] = {'attendees': 0,
                         'user_attending': 0,
                         'user_waitlisted': 0,
                         'max_attendees': row[2],
                         'waitlist_size': row[3]}

    return events
    
if __name__ == '__main__':
    main()
