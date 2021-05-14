from flask import request
from flask_restful import Resource
import scanpy as sc
import os
import sys
import geardb

def tsne_or_umap_present(ana):
  """Return True if tSNE or UMAP plot was calculated for the given analysis."""
  return ana.tsne['tsne_calculated'] == 1 or ana.tsne['umap_calculated'] == 1

class Analyses(Resource):
    """Resource for retrieving all public and private analysis."""

    def get(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)
        dataset = geardb.get_dataset_by_id(dataset_id)

        acollection = geardb.AnalysisCollection()

        acollection.get_all_by_dataset_id(
            user_id=user.id,
            session_id=session_id,
            dataset_id=dataset_id)


        public_tsne = list(
            filter(
              tsne_or_umap_present,
              acollection.public
            )
          )

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
