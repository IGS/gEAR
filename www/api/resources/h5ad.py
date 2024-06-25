from flask import request
from flask_restful import Resource
import os
import geardb

class H5ad(Resource):
    """H5ad Container

    Returns
    -------
    dict
        Column names in .obs
        TODO: Return all relevant h5ad data
    """
    def get(self, dataset_id):
        args = request.args
        analysis_id = args.get('analysis_id')
        session_id = request.cookies.get('gear_session_id')

        import scanpy as sc

        if analysis_id:
            # session_id = request.cookies.get('gear_session_id')
            user = geardb.get_user_from_session_id(session_id)

            ana = geardb.Analysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user.id)
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
        columns = adata.obs.columns.tolist()

        # Do we have these coordinates from the analysis workbench?
        # If so, add these as potential dropdown choices.
        if hasattr(adata, 'obsm'):
            if 'X_tsne' in adata.obsm:
                columns.append('X_tsne_1')
                columns.append('X_tsne_2')
            if 'X_umap' in adata.obsm:
                columns.append('X_umap_1')
                columns.append('X_umap_2')
            if 'X_pca' in adata.obsm:
                columns.append('X_pca_1')
                columns.append('X_pca_2')

        replicate_fields = ["replicate", "BioRep", "TechRep"]
        has_replicates = 0
        for rep in replicate_fields:
            if rep in columns:
                columns.remove(rep)
                has_replicates = 1
        if 'time_point_order' in columns:
            columns.remove('time_point_order')
        # get a map of all levels for each column

        levels = {}
        for col in columns:
            try:
                levels[col] = adata.obs[col].cat.categories.tolist()

                # if there is missing data, add that as a level
                if adata.obs[col].isnull().sum() > 0 and "NA" not in levels[col]:
                    levels[col].append("NA")

                # Drop level if it has more than 50 unique values (i.e. barcodes)
                # Most likely people won't want to filter/sort/plot with them
                if len(levels[col]) > 50:
                    del levels[col]

            except:
                pass
                # If levels are not categorical I don't believe
                # we need to return levels

        return {
            "success": 1,
            "num_obs": adata.n_obs,
            "obs_columns": columns,
            "obs_levels": levels,
            "has_replicates": has_replicates
        }
