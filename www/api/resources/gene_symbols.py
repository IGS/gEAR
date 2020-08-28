from flask import request
from flask_restful import Resource
import scanpy as sc
import os
import geardb

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

        if analysis_id:
            ana = geardb.Analysis(
                id=analysis_id,
                dataset_id=dataset_id,
                session_id=session_id,
                user_id=user.id
            )

            ana.discover_type()
            adata = sc.read_h5ad(ana.dataset_path())
        else:
            ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
            h5_path = ds.get_file_path()
            # Let's not fail if the file isn't there
            if not os.path.exists(h5_path):
                return {
                    "success": -1,
                    "message": "No h5 file found for this dataset"
                }

            adata = sc.read_h5ad(h5_path)

        return {
            "success": 1,
            "gene_symbols": adata.var.gene_symbol.tolist()
        }