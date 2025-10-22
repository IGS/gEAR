import os

import geardb
from flask import request
from flask_restful import Resource

from .common import get_adata_shadow, get_spatial_adata


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

        ds = geardb.get_dataset_by_id(dataset_id)
        if not ds:
            return {
                "success": -1,
                'message': "No dataset found with that ID"
            }

        try:
            if ds.dtype == "spatial":
                adata = get_spatial_adata(analysis_id, dataset_id, session_id)
            else:
                adata = get_adata_shadow(analysis_id, dataset_id, session_id)

        except FileNotFoundError:
            return {
                "success": -1,
                'message': "No file found for this dataset"
            }
        except Exception as e:
            return {
                "success": -1,
                'message': str(e)
            }

        return {
            "success": 1,
            # Convert from categorical to string, and fill in NaN values
            "gene_symbols": adata.var.gene_symbol.astype(str).fillna("None").tolist()
        }
