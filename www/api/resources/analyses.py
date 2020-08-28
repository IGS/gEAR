from flask import request
from flask_restful import Resource
import scanpy as sc
import os
import sys
import geardb


class Analyses(Resource):
    """Resource for retrieving all public and private analysis."""

    def get(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)
        dataset = geardb.get_dataset_by_id(dataset_id)

        # This is commented right now until we work out modeling URL-shared datasets/profiles
        #if not dataset.is_public:
        #   if user.id != dataset.owner_id:
        #        return {
        #            "success": -1,
        #            "message": 'Only the owner can access this dataset.'
        #        }

        acollection = geardb.AnalysisCollection()
        acollection.get_all_by_dataset_id(
            user_id=user.id,
            session_id=session_id,
            dataset_id=dataset_id)

        public_tsne = list(
            filter(
              lambda ana: ana.tsne['tsne_calculated'] == 1,
              acollection.public
            )
          )

        private_tsne = list(
            filter(
              lambda ana: ana.tsne['tsne_calculated'] == 1,
              acollection.user_saved
            )
          )

        return {
            "success": 1,
            "public": public_tsne,
            "private": private_tsne,
        }
