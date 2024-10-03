from flask import request
from flask_restful import Resource
import os
import geardb

from .common import get_adata_shadow

class GeneSymbols(Resource):
    """Gene Symbols

    This resource can share either all genes in a dataset,
    or if a analysis query is passed, all genes in the h5ad
    within that analysis.

    Returns
    -------
    list
        All gene symbols in dataset
    """
    def get(self, dataset_id):
        analysis_id = request.args.get('analysis')
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)

        ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5_path = ds.get_file_path()

        try:
            adata = get_adata_shadow(analysis_id, dataset_id, session_id, h5_path)
        except FileNotFoundError:
            return {
                "success": -1,
                'message': "No h5 file found for this dataset"
            }

        return {
            "success": 1,
            # Convert from categorical to string, and fill in NaN values
            "gene_symbols": adata.var.gene_symbol.astype(str).fillna("None").tolist()
        }
