#!/opt/bin/python3

"""
Used by projection.html, this script gets a list of the available analysis pattern files.
Returns as a list
"""

import cgi, json, os, sys
from pathlib import Path

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

import geardb

PATTERN_BASE_DIR = "/var/www/patterns"

def main():
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')

    if session_id is not None:
        user = geardb.get_user_from_session_id(session_id)
    else:
        user = None

    # Initialize initial select option
    result = [{"id":"", "title":"Choose set:"}]

    base_dir = Path(PATTERN_BASE_DIR)

    for p_file in base_dir.iterdir():
        name = p_file.name    # get basname
        if not name.endswith(".tab"):   # skip non-tab files, like .gitignore
            continue
        title = name.replace('.tab', '').replace('ROWmeta_DIMRED_', '')
        result.append({"id":name, "title":title})

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
