import os
import requests

from flask import request
from flask_restful import Resource

from .common import ANNOTATION_BEDDB_UID, HIGLASS_URL

# TODO: Need to figure out how to do this, as we have no way to search a gene from a tile/tileset.


def get_gene_coords_from_higlass(tileset_id: str, gene_name: str):
    """
    Queries the HiGlass server's text search/tiles endpoint for a specific gene.
    """
    # HiGlass endpoint for fetching specific tile data or searching tags
    # We query the 'tiles' endpoint for the text-index tile (usually tile 0.0 for search)
    query_url = f"{HIGLASS_URL}/suggest/?d={tileset_id}&ac={gene_name}"

    try:
        response = requests.get(query_url)
        response.raise_for_status()
        data: list = response.json()

        if not len(data):
            print(f"No data returned from HiGlass for tileset {tileset_id} and gene {gene_name}.")
            return None

        # list of dicts.  For now just return the first result in case of fuzzy searches (best score)
        # Fields - chr, txStart, txEnd, score, geneName
        return data[0]
    except requests.RequestException as e:
        print(f"Error querying HiGlass for tileset {tileset_id} and gene {gene_name}: {e}")
        return None


class HiGlassGene(Resource):

    def get(self, gene_symbol):

        # Get the assembly from the query parameters
        assembly = request.args.get('assembly', "")

        beddb_uid = ANNOTATION_BEDDB_UID.get(assembly, "")
        if not beddb_uid:
            raise ValueError(
                f"Assembly {assembly} is not supported or does not have a BEDdb UID."
            )

        return get_gene_coords_from_higlass(beddb_uid, gene_symbol)
