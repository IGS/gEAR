#!/opt/bin/python3

"""
Returns the principal investigator (PI) names already given for datasets, so the uploader's
metadata form can suggest them and each PI is spelled one way.

Names come from public datasets and from the user's own datasets. Each name is returned with the
institute most often given with it (from the dataset's contact institute), to tell apart PIs who
share a name.

Optional parameter: session_id (without it, only public datasets are used)

Returns: {"success": 1, "pi_names": [{"name": "Jane Doe", "institute": "IGS"}, ...]}
"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb


def main():
    print('Content-Type: application/json\n\n')
    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')

    user = geardb.get_user_from_session_id(session_id) if session_id else None
    owner_id = user.id if user else -1

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    # GROUP BY ignores case and accents (column collation), so such variants of a name come back as one
    qry = """
        SELECT pi_name, contact_institute, COUNT(*) AS dataset_count
          FROM dataset
         WHERE pi_name IS NOT NULL AND pi_name != ''
           AND marked_for_removal = 0
           AND (is_public = 1 OR owner_id = %s)
         GROUP BY pi_name, contact_institute
         ORDER BY pi_name, dataset_count DESC
    """
    cursor.execute(qry, (owner_id,))

    # Rows for a name are ordered by count, so the first institute seen is the most common one
    pi_names = {}
    for pi_name, institute, _count in cursor:
        pi_names.setdefault(pi_name.lower(), {'name': pi_name, 'institute': institute or ''})

    cursor.close()
    cnx.close()

    print(json.dumps({'success': 1, 'pi_names': list(pi_names.values())}))


if __name__ == '__main__':
    main()
