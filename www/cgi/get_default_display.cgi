#!/opt/bin/python3

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    user_id = form.getvalue('user_id')
    dataset_id = form.getvalue('dataset_id')

    default_display_id = geardb.get_default_display(
      user_id=user_id, dataset_id=dataset_id)

    if default_display_id is None:
      # User owner's default
      dataset = geardb.get_dataset_by_id(id=dataset_id)
      default_display_id = geardb.get_default_display(
        user_id=dataset.owner_id, dataset_id=dataset_id)

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(dict(default_display_id=default_display_id)))

if __name__ == '__main__':
    main()
