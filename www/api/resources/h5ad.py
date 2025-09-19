import geardb
from flask import request
from flask_restful import Resource

from .common import get_adata_shadow, get_spatial_adata


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

        ds = geardb.get_dataset_by_id(dataset_id)
        if not ds:
            return {
                "success": -1,
                'message': "No dataset found with that ID"
            }
        is_spatial = ds.dtype == "spatial"


        try:
            if is_spatial:
                adata = get_spatial_adata(analysis_id, dataset_id, session_id, include_images=False)
            else:
                adata = get_adata_shadow(analysis_id, dataset_id, session_id)
        except FileNotFoundError:
            return {
                "success": -1,
                'message': "No dataset file found."
            }
        except Exception as e:
            return {
                "success": -1,
                'message': str(e)
            }

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

            except Exception as e:
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
