from flask import request
from flask_restful import Resource

import scanpy as sc
import matplotlib.pyplot as plt

import json
import sys
import geardb
import base64
import io

sc.settings.set_figure_params(dpi=100)

class TSNEData(Resource):
    """Resource for retrieving tsne data from an analysis.

    Returns
    -------
      Byte stream image data
    """
    # This endpoint would be a get request to
    # '/api/plot/<SOME_DATASET_ID>/tsne?gene=<SOME_GENE>&analysis=<SOME_ANALYSIS_ID>
    def get(self, dataset_id):
        #print("DEBUG: TSNEData.get() called", file=sys.stderr)
        gene_symbol = request.args.get('gene')
        analysis_id = request.args.get('analysis')
        colorize_by = request.args.get('colorize_by')
        colors = request.args.get('colors')
        x_axis = request.args.get('x_axis')
        y_axis = request.args.get('y_axis')
        session_id = request.cookies.get('gear_session_id')
        user = geardb.get_user_from_session_id(session_id)
        analysis_owner_id = request.args.get('analysis_owner_id')
        sc.settings.figdir = '/tmp/'

        if not gene_symbol or not dataset_id:
            #print("DEBUG: Request needs both dataset id and gene symbol.", file=sys.stderr)
            return {
                "success": -1,
                "message": "Request needs both dataset id and gene symbol."
            }

        dataset = geardb.get_dataset_by_id(dataset_id)

        # If the dataset is not public, make sure the user
        # requesting resource owns the dataset
        # This is commented right now until we work out modeling URL-shared datasets/profiles
        #if not dataset.is_public:
        #    if user.id != dataset.owner_id:
        #        return {
        #            "success": -1,
        #            "message": 'Only the owner can access this dataset.'
        #        }

        if analysis_id:
            # need analysis_type here, but can discover it
            ana = geardb.Analysis(id=analysis_id, dataset_id=dataset_id,
                                  session_id=session_id,
                                  user_id=analysis_owner_id)
            ana.discover_type()
        else:
            ana = geardb.Analysis(type='primary', dataset_id=dataset_id)

        adata = ana.get_adata(backed=True)

        gene_symbols = (gene_symbol,)
        if 'gene_symbol' in adata.var.columns:
            gene_filter = adata.var.gene_symbol.isin(gene_symbols)
            if not gene_filter.any():
                #print("DEBUG: Gene not found.", file=sys.stderr)
                return {
                    'success': -1,
                    'message': 'Gene not found',
                }

            if gene_filter.sum() > 1:
                #print("DEBUG: Multiple Ensembl IDs matched the gene", file=sys.stderr)
                return {
                    "success": -1,
                    "message": f"Multiple Ensembl IDs matched the gene: {gene_symbols[0]}"
                }
        else:
            #print("DEBUG: Missing gene_symbol in adata.var", file=sys.stderr)
            return {
                'success': -1,
                'message': 'Missing gene_symbol in adata.var'
            }


        if analysis_id is None or analysis_id == 'null' or analysis_id == 'undefined':
            #print("DEBUG: working with primary analysis", file=sys.stderr)
            # working on primary data uploaded, so find tSNE_1 and tSNE_2 in obs and build X_tsne
            adata.obsm['X_tsne'] = adata.obs[[x_axis, y_axis]].values

        # We also need to change the adata's Raw var dataframe
        # We can't explicitly reset its index so we reinitialize it with
        # the newer adata object.
        # https://github.com/theislab/anndata/blob/master/anndata/base.py#L1020-L1022
        if adata.raw is not None:
            adata.raw = adata

        # If colorize_by is passed we need to generate that image first, before the index is reset
        #  for gene symbols, then merge them.
        if colorize_by is not None and colorize_by != 'null':
            # were custom colors passed?  the color index is the 'colorize_by' label but with '_colors' appended
            color_idx_name = "{0}_colors".format(colorize_by)

            ## why 2?  Handles the cases of a stringified "{}" or actual keyed JSON
            if colors is not None and len(colors) > 2:
                colors = json.loads(colors)
                adata.uns[color_idx_name] = []

                for idx in adata.obs[colorize_by].cat.categories:
                    adata.uns[color_idx_name].append(colors[idx])

            # the figsize options here (paired with dpi spec above) dramatically affect the definition of the image
            io_fig = plt.figure(figsize=(13, 4))
            spec = io_fig.add_gridspec(ncols=2, nrows=1, width_ratios=[1.1, 1])

            f1 = io_fig.add_subplot(spec[0,0])
            f2 = io_fig.add_subplot(spec[0,1])

            adata.var = adata.var.reset_index().set_index('gene_symbol')
            sc.pl.tsne(adata, color=[gene_symbol], color_map='YlOrRd', ax=f1, show=False)
            sc.pl.tsne(adata, color=[colorize_by], ax=f2, show=False)
            f2.legend(bbox_to_anchor=[1, 1], ncol=2)
        else:
            adata.var = adata.var.reset_index().set_index('gene_symbol')
            io_fig = sc.pl.tsne(adata, color=[gene_symbol], color_map='YlOrRd', return_fig=True)

        io_pic = io.BytesIO()
        io_fig.tight_layout()
        io_fig.savefig(io_pic, format='png')
        io_pic.seek(0)

        return base64.b64encode(io_pic.read()).decode("utf-8")







