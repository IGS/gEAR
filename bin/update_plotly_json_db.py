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
    SELECT * from dataset_display;
"""
read_cursor.execute(get_all_displays)
for (id_, dataset_id, user_id, label, plot_type, plotly_config) in read_cursor:
    if plot_type == 'bar' or plot_type == 'line' or plot_type == 'violin':
      config = json.loads(plotly_config)
      if type(config) is str:
        # This shouldn't happen...but some configs have been double
        # encoded as json.
        config = json.loads(plotly_config)

      index = config.get('index')
      if index and len(index) == 3:
        del config['index']

        facet_row, color_name, x = index
        config['facet_row'] = facet_row if plot_type == "line" else None
        config['facet_col'] = facet_row if plot_type != "line" else None
        config['color_name'] = color_name
        config['x_axis'] = x

      elif index and len(index) == 2:
        del config['index']

        color_name, x = index
        config['color_name'] = color_name
        config['x_axis'] = x,
        config['facet_row'] = None
        config['facet_col'] = None

      elif index and len(index) == 1:
        del config['index']

        x, = index
        config['x_axis'] = x
        config['color_name'] = None
        config['facet_row'] = None
        config['facet_col'] = None

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
