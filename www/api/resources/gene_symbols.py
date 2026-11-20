"""
gene_symbols.py - List the gene symbols in a dataset or analysis.

Serves /h5ad/<dataset_id>/genes in www/api/api.py.
"""

import geardb
from flask import Response, request
from flask_restful import Resource

from .common import get_adata_shadow, get_spatial_adata


def text_response(body: str, status: int = 200) -> Response:
    """Plain-text response shown inline in the browser (e.g. a gene list opened in a new tab)."""
    return Response(
        body,
        status=status,
        mimetype="text/plain",
        headers={"Content-Disposition": "inline"},
    )


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
        """
        Return all gene symbols in the dataset, or in the analysis given by the "analysis_id" query param.

        Query params:
            analysis_id: analysis ID (optional). "analysis" is also accepted, for older callers.
            format: "json" (default) or "text". JSON lists every var row in file order, with missing
                symbols as "None". Text lists one unique symbol per line, sorted case-insensitively,
                without missing symbols.
        """
        analysis_id = request.args.get('analysis_id') or request.args.get('analysis')
        as_text = request.args.get('format', 'json').lower() == 'text'
        session_id = request.cookies.get('gear_session_id')

        ds = geardb.get_dataset_by_id(dataset_id)
        if not ds:
            message = "No dataset found with that ID"
            if as_text:
                return text_response(message, 404)
            return {
                "success": -1,
                'message': message
            }

        try:
            if ds.dtype == "spatial":
                adata = get_spatial_adata(analysis_id, dataset_id, session_id)
            else:
                adata = get_adata_shadow(analysis_id, dataset_id, session_id)

        except FileNotFoundError:
            message = "No file found for this dataset"
            if as_text:
                return text_response(message, 404)
            return {
                "success": -1,
                'message': message
            }
        except Exception as e:
            if as_text:
                return text_response(str(e), 400)
            return {
                "success": -1,
                'message': str(e)
            }

        if as_text:
            symbols = {str(symbol) for symbol in adata.var.gene_symbol.dropna()}
            symbols.discard("None")
            return text_response("\n".join(sorted(symbols, key=lambda s: (s.lower(), s))) + "\n")

        return {
            "success": 1,
            # Fill in NaN values, then convert from categorical to string
            "gene_symbols": adata.var.gene_symbol.astype(object).fillna("None").astype(str).tolist()
        }
