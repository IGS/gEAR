#!/opt/bin/python

# This is to fix an issue where some layouts have duplicated display members, if the user
# saved layouts in the layout arranger when the duplcation bug was active (https://github.com/IGS/gEAR/issues/768)

import sys

from pathlib import Path
lib_path = Path(__file__).resolve().parents[1].joinpath('lib')

sys.path.append(str(lib_path))

import geardb

conn = geardb.Connection()
cursor = conn.get_cursor()

# print row count
qry = "SELECT COUNT(*) FROM layout_displays"
cursor.execute(qry)
row_count = cursor.fetchone()[0]
print("Row count before deletion: {}".format(row_count))

# https://www.tutorialspoint.com/mysql/mysql-delete-duplicate-records.htm
qry = """
DELETE ld1 FROM layout_displays ld1
INNER JOIN layout_displays ld2
WHERE ld1.layout_id = ld2.layout_id
AND ld1.display_id = ld2.display_id
AND ld1.start_col = ld2.start_col
AND ld1.grid_width = ld2.grid_width
AND ld1.start_row = ld2.start_row
AND ld1.grid_height = ld2.grid_height
AND ld1.id > ld2.id
"""
cursor.execute(qry)

cursor.commit()

# print row count
qry = "SELECT COUNT(*) FROM layout_displays"
cursor.execute(qry)
row_count = cursor.fetchone()[0]
print("Row count after deletion: {}".format(row_count))

cursor.close()
conn.close()