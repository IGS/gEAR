import os

import geardb
from flask import request
from flask_restful import Resource
from gear.analysis import get_analysis

from .common import get_adata_from_analysis, get_spatial_adata


class TopPCAGenes(Resource):
    """Plot Top Genes of Prinicpal Components"""
    def post(self):
        req = request.form
        analysis_id = req.get('analysis_id')
        analysis_type = req.get('analysis_type')
        pcs = [int(pc.replace(' ', '')) for pc in req.get('pcs').split(',')]
        # Currently, ScanPy is throwing errors both when number
        # of PCs is less than 2 or greater than 5. Have an issue
        # open: https://github.com/theislab/scanpy/issues/431
        if len(pcs) < 2 or len(pcs) > 5:
            return {
                "success": 0,
                "message": "Number of PCs must between 2 and 5",
            }

        dataset_id = req.get('dataset_id')
        session_id = req.get('session_id')

        if not dataset_id:
            return {
                "success": 0,
                "message": "Missing dataset_id",
            }

        ds = geardb.get_dataset_by_id(dataset_id)
        if not ds:
            return {
                "success": 0,
                "message": "Invalid dataset_id",
            }

        is_spatial = ds.dtype == "spatial"

        analysis = None
        if analysis_id or analysis_type:
            analysis = {}
            if analysis_id:
                analysis['id'] = analysis_id
            if analysis_type:
                analysis['type'] = analysis_type

        ana = get_analysis(analysis, dataset_id, session_id, is_spatial=is_spatial)

        dest_datafile_path = ana.dataset_path
        dest_directory = os.path.dirname(dest_datafile_path)

        #print("DEBUG: dest_directory: {0}".format(dest_directory), file=sys.stderr)

        if not os.path.exists(dest_directory):
            os.makedirs(dest_directory)

        try:
            adata = ana.get_adata()
        except Exception as e:
            return {
                "success": -1,
                "message": str(e),
            }


        import scanpy as sc
        sc.settings.figdir = dest_directory + "/figures"

        #print("DEBUG: Writing file to path: {0}".format(dest_datafile_path), file=sys.stderr)

        adata.var = adata.var.reset_index().set_index('gene_symbol')

        try:
            sc.pl.pca_loadings(adata, components=pcs, save='.png')
        except Exception as e:
            return {
                "message": str(e),
                "success": -1,
            }
        # Code below is to grab top genes from PCA. Will have to use
        # this in order to get a heatmap implementation working.
        # sort_idcs = np.argsort(adata.varm['PCs'][:, 0])
        # genes_ranked_by_loading_in_PC1 = adata.var_names[sort_idcs]
        return {
            "success": 1,
        }
