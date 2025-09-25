#!/opt/bin/python3

"""
get_PCs_from_anndata.cgi

Description: Given the current AnnData object
1. Retrieve the PCs stored in the "varm" area (numpy array)
2. Add the genes as the row indexes and PC names as column labels
3. Return as JSON
"""

import cgi
import json
import os
import sys
from pathlib import Path

import pandas as pd

original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

TWO_LEVELS_UP = 2
abs_path_gear = Path(__file__).resolve().parents[TWO_LEVELS_UP]
abs_path_lib = abs_path_gear.joinpath('lib')
# abs_path_lib is a Path object so we need to convert to string
sys.path.insert(0, str(abs_path_lib))

import geardb
from gear.analysis import get_analysis


def main():
    form = cgi.FieldStorage()
    analysis_id = form.getvalue('analysis_id')
    analysis_type = form.getvalue('analysis_type')
    dataset_id = form.getvalue('dataset_id')
    session_id = form.getvalue('session_id')

    if not dataset_id:
        return_error_response("No dataset_id provided")
        return

    ds = geardb.get_dataset_by_id(dataset_id)
    if not ds:
        return_error_response("No dataset found with that ID")
        return

    is_spatial = ds.dtype == "spatial"

    analysis_obj = dict(id=analysis_id, type=analysis_type) if analysis_id and analysis_type else None

    try:
        ana = get_analysis(analysis_obj, dataset_id, session_id, is_spatial=is_spatial)
    except Exception:
        return_error_response("Analysis for this dataset is unavailable.")

    try:
            args = {}
            if is_spatial:
                args['include_images'] = False
            adata = ana.get_adata(**args)
    except Exception:
        return_error_response("Could not create dataset object using analysis.")

    if "PCs" not in adata.varm:
        return_error_response("PCs not found in AnnData object")
        return

    pc_data = adata.varm["PCs"]
    columns = make_pc_columns(pc_data.shape[1])
    genes = adata.var.index.tolist()
    gene_symbols = adata.var.gene_symbol.tolist()    # Not needed for the weighted gene cart but are needed to make the geardb.Gene objects

    df = pd.DataFrame(pc_data, index=genes, columns=columns) # type: ignore

    json_obj = json.loads(df.to_json(orient="split"))

    result = {'success': 1, "message": "Retrieved PCs for genes", "pc_data":json_obj, "gene_symbols":gene_symbols}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

def make_pc_columns(num_pcs):
    columns = []
    for i in range(num_pcs):
        columns.append("PC" + str(i + 1))
    return columns

def return_error_response(msg):
    result = dict()
    result['success'] = 0
    result['message'] = msg
    result["pc_data"] = None
    result["gene_symbols"] = None
    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()