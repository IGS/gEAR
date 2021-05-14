from flask import request
from flask_restful import Resource
import scanpy as sc
import os, sys
import geardb

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

        user = geardb.get_user_from_session_id(session_id)
        ana = geardb.Analysis(
            id=analysis_id,
            type=analysis_type,
            dataset_id=dataset_id,
            session_id=session_id,
            user_id=user.id
        )

        source_datafile_path = ana.dataset_path()

        dest_datafile_path = ana.dataset_path()
        dest_directory = os.path.dirname(dest_datafile_path)

        #print("DEBUG: dest_directory: {0}".format(dest_directory), file=sys.stderr)

        if not os.path.exists(dest_directory):
            os.makedirs(dest_directory)

        adata = sc.read_h5ad(source_datafile_path)
        sc.settings.figdir = dest_directory + "/figures"
        
        #print("DEBUG: Writing file to path: {0}".format(dest_datafile_path), file=sys.stderr)

        adata.var = adata.var.reset_index().set_index('gene_symbol')

        try:
            sc.pl.pca_loadings(adata, components=pcs, save='.png')
        except:
            return {
                "success": -1,
            }
        # Code below is to grab top genes from PCA. Will have to use
        # this in order to get a heatmap implementation working.
        # sort_idcs = np.argsort(adata.varm['PCs'][:, 0])
        # genes_ranked_by_loading_in_PC1 = adata.var_names[sort_idcs]
        return {
            "success": 1,
        }
