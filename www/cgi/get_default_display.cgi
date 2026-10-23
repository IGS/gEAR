#!/opt/bin/python3

"""
get_default_display.cgi - Get the default display ID for a dataset.

Input: dataset_id, session_id (optional), is_multigene (0/1).
Output: JSON {default_display_id, grid_width, width_rule, width_counts}; falls back to the dataset
        owner's default when the user has none set. grid_width is the suggested width (of 12) for
        automatically arranged layouts, or null if there is no basis for one. width_rule says why
        (see geardb.get_auto_layout_width) and width_counts is how many collections used each width.
"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    dataset_id = form.getfirst('dataset_id')
    is_multigene = int(form.getfirst('is_multigene', 0))

    user = geardb.get_user_from_session_id(session_id=session_id)

    default_display_id = None
    if user:
      default_display_id = geardb.get_default_display(
        user_id=user.id, dataset_id=dataset_id, is_multigene=is_multigene
)
    dataset = None
    if default_display_id is None:
      dataset = geardb.get_dataset_by_id(d_id=dataset_id)
    if dataset is not None:
      # User owner's default
      default_display_id = geardb.get_default_display(
        user_id=dataset.owner_id
        , dataset_id=dataset_id
        , is_multigene=is_multigene
      )

    grid_width = None
    width_rule = None
    width_counts = {}
    if default_display_id is not None:
      choice = geardb.get_auto_layout_width(default_display_id)
      grid_width = choice['width']
      width_rule = choice['rule']
      width_counts = {str(width): count for width, count in choice['width_counts'].items()}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(dict(default_display_id=default_display_id, grid_width=grid_width, width_rule=width_rule, width_counts=width_counts)))

if __name__ == '__main__':
    main()
