#!/usr/bin/env python3

import argparse
import mysql.connector
import json

parser = argparse.ArgumentParser( description='')
parser.add_argument('-u', '--user', type=str, required=True)
parser.add_argument('-p', '--password', type=str, required=True)
parser.add_argument('-d', '--database', type=str, required=True)

args = parser.parse_args()
cnx = mysql.connector.connect(user=args.user, password=args.password, database=args.database)

read_cursor = cnx.cursor(buffered=True)
update_cursor = cnx.cursor(buffered=True)

get_all_displays = """
  select * from dataset_display;
"""

read_cursor.execute(get_all_displays)
for (id_, dataset_id, user_id, label, plot_type, plotly_config) in read_cursor:
    config = json.loads(plotly_config)
    if type(config) is str:
      config = json.loads(config)

    update_config_query = """
      UPDATE dataset_display
      SET plotly_config = %s
      WHERE id = %s
    """
    update_cursor.execute(update_config_query, (json.dumps(config), id_))
    cnx.commit()

update_cursor.close()
read_cursor.close()
cnx.close()
