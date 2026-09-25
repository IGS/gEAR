"""
analyses.py - List public and user-saved analyses for a dataset.

Serves /h5ad/<dataset_id>/analyses in www/api/api.py.
"""

from flask import request
from flask_restful import Resource
import geardb
from gear.analysis import AnalysisCollection

def tsne_or_umap_present(ana):
  """Return True if tSNE or UMAP plot was calculated for the given analysis."""
  return ana.tsne['tsne_calculated'] == 1 or ana.tsne['umap_calculated'] == 1

class Analyses(Resource):
    """Resource for retrieving all public and private analysis."""

    def get(self, dataset_id):
        """
        Return public and private analyses for the dataset that have tSNE or UMAP calculated.

        The user is identified by the gear_session_id cookie; "private" is empty when
        not logged in.
        """
        session_id = request.cookies.get('gear_session_id', None)
        user = geardb.get_user_from_session_id(session_id)

        user_id = user.id if user else None

        ds = geardb.get_dataset_by_id(dataset_id)
        if not ds:
            return {
                "success": -1,
                'message': "No dataset found with that ID"
            }

        is_spatial = ds.dtype == "spatial"

        acollection = AnalysisCollection()

        acollection.get_all_by_dataset_id(
            user_id=user_id,
            session_id=session_id,
            dataset_id=dataset_id,
            is_spatial=is_spatial
        )


        public_tsne = list(
            filter(
              tsne_or_umap_present,
              acollection.public
            )
          )

        # This will be empty if the user is not logged in (handled in get_all_by_dataset_id)
        private_tsne = list(
            filter(
              tsne_or_umap_present,
              acollection.user_saved
            )
          )

        return {
            "success": 1,
            "public": public_tsne,
            "private": private_tsne,
        }
