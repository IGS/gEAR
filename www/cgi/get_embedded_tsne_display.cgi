#!/opt/bin/python3

"""
Display entries are usually kept in the database, but for some analysis types (like embedded
tSNE) a database entry isn't necessary, as the dataset itself contains the definition.  What
we need here is a tool which will check for the columns involved and return a JSON structure
as if it were a plotly_config stored in the database.

Needs to return a structure like:

{
  'plotly_config': {
                    'colorize_legend_by': 'cell_type',
                    'x_axis': 'tSNE_1',
                    'y_axis': 'tSNE_2'
                   }
}

"""


import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

# These should match in bin/add_primary_analyses_to_datasets.py
VALID_CLUSTER_COLUMN_NAMES = ['cluster', 'cell_type', 'subclass_label', 'cluster_label', 'subclass_label', 'joint_cluster_round4_annot']
VALID_TSNE_PAIRS = [['tSNE_1', 'tSNE_2'], ['tSNE1', 'tSNE2'], ['tsne1_combined', 'tsne2_combined'],
                    # These are all carlo's custom ones.  Need to resolve this a different way later
                    ['PC1%var6.14', 'PC2%var1.79']
]

def main():
    original_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')
    if not dataset_id:
        sys.stdout = original_stdout
        print('Content-Type: application/json\n\n')
        print(json.dumps({'error': 'No dataset_id provided.'}))
        return

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        print('Content-Type: application/json\n\n')
        print(json.dumps({'error': 'No dataset found with that ID.'}))
        return

    is_spatial = ds.dtype == "spatial"

    analysis_obj = dict(type="primary")

    try:
        ana = geardb.get_analysis(analysis_obj, dataset_id, None, is_spatial=is_spatial)
    except Exception:
        print('Content-Type: application/json\n\n')
        print(json.dumps({'error': 'Could not retrieve analysis'}))
        return

    try:
            args = {}
            if is_spatial:
                args['include_images'] = False
            else:
                args['backed'] = True
            adata = ana.get_adata(**args)
    except Exception:
        print('Content-Type: application/json\n\n')
        print(json.dumps({'error': 'Could not retrieve AnnData object.'}))
        return

    results = {
        'plotly_config': {
            'colorize_legend_by': None,
            'x_axis': None,
            'y_axis': None
        }
    }

    cols = adata.obs.columns.tolist()

    for vname in VALID_CLUSTER_COLUMN_NAMES:
        if vname in cols:
            results['plotly_config']['colorize_legend_by'] = vname # type: ignore
            break

    for pair in VALID_TSNE_PAIRS:
        if pair[0] in cols and pair[1] in cols:
            results['plotly_config']['x_axis'] = pair[0] # type: ignore
            results['plotly_config']['y_axis'] = pair[1] # type: ignore
            break

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(results))

if __name__ == '__main__':
    main()
