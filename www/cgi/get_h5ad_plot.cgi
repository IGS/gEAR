#!/opt/bin/python3

"""

# testing profile:
select id, dtype, organism_id from dataset where id in ('64485ca3-cf99-2993-99a3-54df3a09195c', '6fdd350c-4f82-07e2-3a39-408f105db16d', 'liasdf97-e9a2-po1u-kj11-1k282bjg8j81');
+--------------------------------------+-----------------------+-------------+
| id                                   | dtype                 | organism_id |
+--------------------------------------+-----------------------+-------------+
| 64485ca3-cf99-2993-99a3-54df3a09195c | linegraph-standard    |           5 |
| 6fdd350c-4f82-07e2-3a39-408f105db16d | svg-expression        |           1 |
| liasdf97-e9a2-po1u-kj11-1k282bjg8j81 | image-static-standard |           1 |
+--------------------------------------+-----------------------+-------------+

"""

import cgi
import json
import os
import sys

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDERR, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear.plotting import PlotFactory, PlotDataFactory, get_available_plot_types
import scanpy as sc
sc.settings.verbosity = 0

import pandas as pd

def main():
    result = { 'success': 0 }

    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')
    gene_symbol = form.getvalue('gene_symbol')
    session_id = form.getvalue('session_id')
    group_by = form.getvalue('group_by')
    if group_by:
        group_by = json.loads(group_by)
    colors = form.getvalue('colors')

    if colors:
        colors = json.loads(colors)
    else:
        colors = dict()

    ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
    h5_path = ds.get_file_path()

    # Let's not fail if the file isn't there
    if not os.path.exists(h5_path):
        result['success'] = 0
        result['error'] = "No h5 file found for this dataset"
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        sys.exit()

    adata = sc.read_h5ad(h5_path)
    gene_symbols = (gene_symbol,)

    if 'gene_symbol' in adata.var.columns:

        gene_filter = adata.var.gene_symbol.isin(gene_symbols)
        if not gene_filter.any():
            result['success'] = -1
            result['error'] = 'Gene not found'
            sys.stdout = original_stdout
            print('Content-Type: application/json\n\n')
            print(json.dumps(result))
            sys.exit()
    else:
        # TODO
        # look up ensembl id for gene symbols
        # gene_filter = adata.var.index.isin((ensembl_ids,))
        result['success'] = 0
        result['error'] = "Missing gene_symbol column in adata.var"
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps(result))
        sys.exit()

    selected = adata[:, gene_filter]
    df = pd.DataFrame({"raw_value":selected.X, **selected.obs})

    plot_types_dicts = get_available_plot_types(df, index=group_by)
    result['plot_types'] = plot_types_dicts # ex: [{plot_type: 'violin'}]
    plot_types = list(map(lambda d: d['plot_type'], plot_types_dicts))
    plots = [PlotFactory.create_plot(plot_type, data=df, group_by=group_by, colors=colors) for plot_type in plot_types]
    result['plot_json'] = {plot_type:json.loads(plot.render()) for (plot_type, plot) in zip(plot_types, plots)}
    # plot has to be rendered before colors attribute is updated with default colors
    # TODO: Fix this by generating default colors on initialization
    result['plot_colors'] = plots[0].colors
    result['plot_config'] = PlotFactory.get_config()
    if not group_by:
        group_by = PlotDataFactory.get_groupby_conditions(df)
    result['plot_groups'] = group_by

    result['success'] = 1

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
