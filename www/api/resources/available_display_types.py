from flask import request
from flask_restful import Resource
import os
import geardb

def tsne_or_umap_present(ana):
  """Return True if tSNE or UMAP plot was calculated for the given analysis."""
  return ana.tsne['tsne_calculated'] == 1 or ana.tsne['umap_calculated'] == 1


class MGAvailableDisplayTypes(Resource):
    """Resource for retrieving available multigene display types from given dataset id

    Parameters
    ----------
    dataset_id: str
      Dataset ID
    session_id: str
      Session ID
    analysis_id: str
      Analysis ID

    Returns
    -------
    dict
        Available display types

    """

    def post(self, dataset_id):
        req = request.get_json()
        dataset_id = req.get('dataset_id')
        session_id = req.get('session_id')
        analysis_id = req.get('analysis_id')

        dotplot = True
        heatmap = True
        mg_violin = True
        volcano = True
        quadrant = True

        ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5_path = ds.get_file_path()

        # Import here so that Flask-RESTful does not import it with every API call.
        import scanpy as sc

        # Have a public dataset or user_saved dataset
        if analysis_id:
            user = geardb.get_user_from_session_id(session_id)

            ana = geardb.Analysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user.id)
            ana.discover_type()

            ### Not needed for MG yet

        else:
          # Dataset is primary type
          if not os.path.exists(h5_path):
              # Let's not fail if the file isn't there
              return {
                  "success": -1,
                  'message': "No h5 file found for this dataset"
              }
          adata = sc.read_h5ad(h5_path)

        columns = adata.obs.columns.tolist()

        if "replicate" in columns:
          columns.remove('replicate')

        columns = [col for col in columns if not col.endswith('_colors')]

        # Filter only categorical columns
        categorical_columns = list(filter(lambda x: adata.obs[x].dtype.name == 'category', columns))

        # If there are no categorical columns, then we can't draw any plots so return error
        if len(categorical_columns) == 0:
          return {
            "success": -1,
            'message': "No categorical columns found in this dataset, so plots cannot be drawn. Please choose another dataset."
          }

        # Volcano plots must have at least 2 distinct categorical columns
        if len(categorical_columns) < 2:
          volcano = False

        # Quadrant plots must have at least 3 distinct categorical columns
        if len(categorical_columns) < 3:
          quadrant = False

        available_display_types = {
          "dotplot": dotplot,
          "heatmap": heatmap,
          "mg_violin": mg_violin,
          "volcano": volcano,
          "quadrant": quadrant
        }

        return available_display_types, 200

class AvailableDisplayTypes(Resource):
    """Resource for retrieving available display types from given dataset id

    Parameters
    ----------
    dataset_id: str
      Dataset ID
    session_id: str
      Session ID
    analysis_id: str
      Analysis ID

    Returns
    -------
    dict
        Available display types
    """
    def post(self, dataset_id):
        req = request.get_json()
        dataset_id = req.get('dataset_id')
        session_id = req.get('session_id')
        analysis_id = req.get('analysis_id')

        line = False
        scatter = True  # Should always have scatter as an option
        #contour = False
        tsne_static = False
        umap_static = False
        pca_static = False
        tsne_umap_pca_dynamic = False
        svg_exists = False


        ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5_path = ds.get_file_path()
        # Determine if SVG exists for the primary dataset.
        (base_path, _) = h5_path.split('/datasets/')
        svg_path = f"{base_path}/datasets_uploaded/{dataset_id}.svg"
        # santize svg_path to prevent path traversal
        base_path = os.path.normpath(base_path)
        full_svg_path = os.path.normpath(svg_path)
        if full_svg_path.startswith(base_path) and os.path.exists(full_svg_path):
          svg_exists = True

        # Import here so that Flask-RESTful does not import it with every API call.
        import scanpy as sc

        # Have a public dataset or user_saved dataset
        if analysis_id:
            # session_id = request.cookies.get('gear_session_id')
            user = geardb.get_user_from_session_id(session_id)

            ana = geardb.Analysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user.id)
            ana.discover_type()

            adata = sc.read_h5ad(ana.dataset_path())
            if hasattr(adata, 'obsm') and 'X_tsne' in adata.obsm:
                tsne_static = True
                tsne_umap_pca_dynamic = True
            elif hasattr(adata, 'obsm') and 'X_umap' in adata.obsm:
                umap_static = True
                tsne_umap_pca_dynamic = True
            elif hasattr(adata, 'obsm') and 'X_pca' in adata.obsm:
                pca_static = True
                tsne_umap_pca_dynamic = True
        else:
          # Dataset is primary type
          if not os.path.exists(h5_path):
              # Let's not fail if the file isn't there
              return {
                  "success": -1,
                  'message': "No h5 file found for this dataset"
              }
          adata = sc.read_h5ad(h5_path)

        columns = adata.obs.columns.tolist()

        if len(columns) == 0:
          return {
            "success": -1,
            'message': "No metadata columns found in this dataset, so plots cannot be drawn. Please choose another dataset."
          }

        if "replicate" in columns:
          columns.remove('replicate')

        columns = [col for col in columns if not col.endswith('_colors')]

        for col in columns:
          # If any of the columns are of numerical type,
          # then we can draw line plots
          if "float" in str(adata.obs[col].dtype) or "int" in str(adata.obs[col].dtype):
            line = True

          """
          # if any columns have tSNE in the name, enable support for that.
          if 'tsne'.lower() in str(col).lower():
            tsne_static = True
            tsne_umap_pca_dynamic = True

          # if any columns have UMAP in the name, enable support for that.
          if 'umap'.lower() in str(col).lower():
            umap_static = True
            tsne_umap_pca_dynamic = True

          # if any columns have PCA in the name, enable support for that.
          if 'pca'.lower() in str(col).lower():
            pca_static = True
            tsne_umap_pca_dynamic = True

          # Carlo wants to be able to plot DimRed columns this way
          if str(col).startswith('DimRed') or str(col).startswith('PC'):
            tsne_static = True
            tsne_umap_pca_dynamic = True
          """

        # if at least two columns are float or int, enable tsne/umap/pca plots
        if len([col for col in columns if "float" in str(adata.obs[col].dtype) or "int" in str(adata.obs[col].dtype)]) >= 2:
          tsne_umap_pca_dynamic = True
          tsne_static = True
          umap_static = True
          pca_static = True

        available_display_types = {
          "scatter": scatter,
          "tsne_static": tsne_static,
          "umap_static": umap_static,
          "pca_static": pca_static,
          "tsne/umap_dynamic": tsne_umap_pca_dynamic,
          # Some single cell datasets might have no columns
          "bar": True if len(columns) else False,
          "violin": True if len(columns) else False,
          "line": True if line or "time_point" in columns else False,
          #"contour":True if len(columns) else False,
          "svg": svg_exists
        }

        return available_display_types, 200
