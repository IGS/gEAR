from flask import request
from flask_restful import Resource
import scanpy as sc
import os
import sys
import geardb

class AvailableDisplayTypes(Resource):
    """Resource for retrieving available display types from given dataset id

    Parameters
    ----------
    user_id: str
      User ID
    dataset_id: str
      Dataset ID
    session_id: str
      Session ID

    Returns
    -------
    dict
        Available display types
    """
    def post(self, dataset_id):
        req = request.get_json()
        user_id = req.get('user_id')
        dataset_id = req.get('dataset_id')
        session_id = req.get('session_id')
        analysis_id = req.get('analysis_id')

        line = False
        scatter = False
        tsne_static = False
        tsne_dynamic = False

        if analysis_id:
            # session_id = request.cookies.get('gear_session_id')
            user = geardb.get_user_from_session_id(session_id)

            ana = geardb.Analysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user.id)
            ana.discover_type()
            adata = sc.read_h5ad(ana.dataset_path())
            if hasattr(adata, 'obsm') and hasattr(adata.obsm, 'X_tsne'):
                tsne_static = True
                tsne_dynamic = True

        acollection = geardb.AnalysisCollection()
        acollection.get_all_by_dataset_id(
          user_id=user_id,
          session_id=session_id,
          dataset_id=dataset_id)

        ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5_path = ds.get_file_path()
        # Let's not fail if the file isn't there
        if not os.path.exists(h5_path):
            return {
                "success": -1,
                'message': "No h5 file found for this dataset"
            }

        (base_path, _) = h5_path.split('/datasets/')
        svg_path = f"{base_path}/datasets_uploaded/{dataset_id}.svg"

        # In case these keys don't exist in ana
        try:
          public_tsne = list(
            map(
              lambda ana: ana.tsne['tsne_calculated'],
              acollection.public
            )
          )
        except:
          public_tsne = []

        try:
          private_tsne = list(
            map(
              lambda ana: ana.tsne['tsne_calculated'],
              acollection.user_saved
            )
          )
        except:
          private_tsne = []

        adata = sc.read_h5ad(h5_path)
        columns = adata.obs.columns.tolist()

        if "replicate" in columns:
          columns.remove('replicate')

        for col in columns:
          # If any of the columns are of numerical type,
          # then we can draw both scatter plots and line plots
          if "float" in str(adata.obs[col].dtype) or "int" in str(adata.obs[col].dtype):
            scatter = True
            line = True

          # if any columns have tSNE in the name, enable support for that.
          if 'tsne' in str(col).lower():
              tsne_static = True
              tsne_dynamic = True

          # Carlo wants to be able to plot DimRed columns this way
          if str(col).startswith('DimRed') or str(col).startswith('PC'):
              tsne_static = True
              tsne_dynamic = True

        available_display_types = {
          "scatter": scatter,
          "tsne_static": tsne_static,
          "tsne_dynamic": tsne_dynamic,
          # Some single cell datasets might have no columns
          "bar": True if len(columns) else False,
          "violin": True if len(columns) else False,
          "line": True if line or "time_point" in columns else False,
          "svg": os.path.exists(svg_path)
        }

        # tsne is available if there are any public or private
        # analysis associated with the dataset or the user
        if any(public_tsne) or any(private_tsne):
            available_display_types['tsne_static'] = True
            available_display_types['tsne_dynamic'] = True

        return available_display_types, 200
