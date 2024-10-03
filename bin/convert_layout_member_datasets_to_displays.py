#!/opt/bin/python3

# This script is used to convert layout member dataset IDs into their corresponding
# default single- and multi-gene display IDs.  No arguments are required but it must
# be run from the bin directory.

# This assumes that the dataset-based layout_member table is the original layout member table
# and that a new display-based "layout_member_new" table has been created.

# ! This will not wipe the layout_displays table, so if you run this script multiple times, it will create duplicate records !

import sys

# This resolves some "no localization support" error
from mysql.connector.locales.eng import client_error

from pathlib import Path
lib_path = Path(__file__).resolve().parents[1].joinpath('lib')

sys.path.append(str(lib_path))

import geardb

conn = geardb.Connection()
cursor = conn.get_cursor()

# Get all layout_members

print("Retrieving layout members...", file=sys.stderr)
layout_member_qry = """
              SELECT lm.layout_id, lm.dataset_id, lm.grid_position, lm.mg_grid_position,
              lm.start_col, lm.mg_start_col, lm.grid_width, lm.mg_grid_width,
              lm.start_row, lm.mg_start_row, lm.grid_height, lm.mg_grid_height
              FROM layout_members lm
              JOIN dataset ds ON lm.dataset_id=ds.id
              WHERE ds.marked_for_removal = 0
              ORDER BY lm.layout_id, lm.grid_position
              """


cursor.execute(layout_member_qry)

layout_members = []

for row in cursor:

    layout_id = row[0]
    dataset_id = row[1]
    grid_position = row[2]
    mg_grid_position = row[3]
    start_col = row[4]
    mg_start_col = row[5]
    grid_width = row[6]
    mg_grid_width = row[7]
    start_row = row[8]
    mg_start_row = row[9]
    grid_height = row[10]
    mg_grid_height = row[11]

    if mg_grid_position == 0:
        mg_grid_position = grid_position

    layout_members.append([layout_id, dataset_id, grid_position, mg_grid_position, start_col, mg_start_col, grid_width, mg_grid_width, start_row, mg_start_row, grid_height, mg_grid_height])

# Determine if we are in "legacy" mode where every member start_col is 1
legacy = False
if len(layout_members) > 0:
    # If all members have start_col = 1, then we are in legacy mode
    if all([m[4] == 1 for m in layout_members]):
        #print("Legacy mode found ... rebuilding layout member grid positions...", file=sys.stderr)
        legacy = True

# If "legacy" mode, adjust the start_col and start_row, as well as mg_start_col and mg_start_row
if legacy:
    print("Legacy mode found ... rebuilding layout member grid positions...", file=sys.stderr)
    current_col = 1
    current_row = 1
    curr_layout = -1
    for lm in layout_members:
        #start_col = lm[4]
        #grid_width = lm[6]
        #start_row = lm[8]
        if curr_layout != lm[0]:
            curr_layout = lm[0]
            current_col = 1
            current_row = 1

        width = lm[6]
        if current_col + width > 13:
            current_col = 1
            current_row += 1
        lm[4] = current_col
        lm[8] = current_row
        current_col += width

    # Multigene displays now
    current_col = 1
    current_row = 1
    curr_layout = -1
    for lm in layout_members:
        #mg_start_col = lm[5]
        #mg_grid_width = lm[7]
        #mg_start_row = lm[9]
        if curr_layout != lm[0]:
            curr_layout = lm[0]
            current_col = 1
            current_row = 1

        width = lm[7]
        if current_col + width > 13:
            current_col = 1
            current_row += 1
        lm[5] = current_col
        lm[9] = current_row
        current_col += width


print("Converted layout members to layout displays...", file=sys.stderr)
for lm in layout_members:

    layout_id = lm[0]
    dataset_id = lm[1]
    grid_position = lm[2]
    mg_grid_position = lm[3]
    start_col = lm[4]
    mg_start_col = lm[5]
    grid_width = lm[6]
    mg_grid_width = lm[7]
    start_row = lm[8]
    mg_start_row = lm[9]
    grid_height = lm[10]
    mg_grid_height = lm[11]


    # get the default single gene and multi gene display ID preferences for the dataset ID

    single_gene_preference = """
        SELECT display_id FROM dataset_preference
        WHERE dataset_id = %s AND is_multigene = 0
        AND user_id = (SELECT owner_id FROM dataset WHERE id = %s)
    """
    multi_gene_preference = """
        SELECT display_id FROM dataset_preference
        WHERE dataset_id = %s AND is_multigene = 1
        AND user_id = (SELECT owner_id FROM dataset WHERE id = %s)
    """

    # Caveats with this approach:
    # 1. If the dataset has no default single- or multi-gene display preference, the layout member will not be inserted
    # 2. Only display preferences by the dataset owner will be used

    cursor.execute(single_gene_preference, (dataset_id,dataset_id))
    single_fetch = cursor.fetchone()

    cursor.execute(multi_gene_preference, (dataset_id,dataset_id))
    multi_fetch = cursor.fetchone()

    # insert the new layout member record as separate single- and multi-gene display records
    new_layout_member_qry = """
        INSERT INTO layout_displays (layout_id, display_id, grid_position, start_col, grid_width, start_row, grid_height, math_preference)
        VALUES (%s, %s, %s, %s, %s, %s, %s, 'raw')
    """
    if single_fetch:
        single_gene_display_id = single_fetch[0]
        cursor.execute(new_layout_member_qry, (layout_id, single_gene_display_id, grid_position, start_col, grid_width, start_row, grid_height))

    if multi_fetch:
        multi_gene_display_id = multi_fetch[0]
        cursor.execute(new_layout_member_qry, (layout_id, multi_gene_display_id, mg_grid_position, mg_start_col, mg_grid_width, mg_start_row, mg_grid_height))

cursor.close()

# commit the changes
conn.commit()
conn.close()
