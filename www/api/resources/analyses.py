from flask import request
from flask_restful import Resource
import geardb

def tsne_or_umap_present(ana):
  """Return True if tSNE or UMAP plot was calculated for the given analysis."""
  return ana.tsne['tsne_calculated'] == 1 or ana.tsne['umap_calculated'] == 1

class Analyses(Resource):
    """Resource for retrieving all public and private analysis."""

    def get(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)

        user_id = user.id if user else None

        acollection = geardb.AnalysisCollection()

        acollection.get_all_by_dataset_id(
            user_id=user_id,
            session_id=session_id,
            dataset_id=dataset_id)


        public_tsne = list(
            filter(
              tsne_or_umap_present,
              acollection.public
            )
          )

        private_tsne = []
        if user_id:
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
