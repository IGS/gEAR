#!/usr/bin/env python3

"""
Some plotly dataset displays are missing important attributes in their
config json. This script will compute what used previously as the index
and will add the attributes to their config.
"""

import argparse
import scanpy as sc
import mysql.connector
import json
import os

# curr_dir_path = os.path.dirname(os.path.abspath(__file__))

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
        config = json.loads(config)


      x_axis = config.get('x_axis')
      color_name = config.get('color_name')
      facet_row = config.get('facet_row')
      facet_col = config.get('facet_col')

      if all([not x_axis, not color_name, not facet_row, not facet_col]):
        if 'index' in config:
          del config['index']
        # If all of these are True, then none of them exist in a user's config
        # assuming we are in the bin directory
        file_path = f"../www/datasets/{dataset_id}.h5ad"
        adata = sc.read(file_path)
        if not adata.obs.empty:

          excluded_columns = ['raw_value', 'replicate', 'time_unit', 'time_point_order']
          filter_cols = ~adata.obs.columns.isin(excluded_columns)
          conditions = (
              adata.obs.loc[:, filter_cols]
                  .nunique()
                  .sort_values()
                  .index
                  .tolist()
          )
          top_conditions = conditions[-3:]
          if len(top_conditions) == 1:
              x_axis, = top_conditions
              config['x_axis'] = x_axis
              config['color_name'] = None
              config['facet_row'] = None
              config['facet_col'] = None
          elif len(top_conditions) == 2:
              color_name, x_axis = top_conditions
              config['color_name'] = color_name
              config['x_axis'] = x_axis
              config['facet_row'] = None
              config['facet_col'] = None
          elif len(top_conditions) == 3:
              facet_row, color_name, x_axis = top_conditions
              config['facet_row'] = facet_row if plot_type == "line" else None
              config['facet_col'] = facet_row if plot_type != "line" else None
              config['color_name'] = color_name
              config['x_axis'] = x_axis
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
