import os
import sys

import geardb
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from flask import request
from flask_restful import Resource
from gear.plotting import PlotError

from .common import clip_expression_values, create_projection_adata


class SvgData(Resource):
    """Resource for retrieving data from h5ad to be used to color svgs.

    Returns
    -------
    dict
        SVG data
    """
    def get(self, dataset_id):
        gene_symbol = request.args.get('gene', None)
        projection_id = request.args.get('projection_id', None)    # projection id of csv output
        expression_min_clip = request.args.get('expression_min_clip', None)  # minimum expression value to clip to, if applicable
        vmin = request.args.get('vmin', None)
        vmax = request.args.get('vmax', None)

        add_user_defined = False
        if vmax is not None or vmin is not None:
            add_user_defined = True

        if expression_min_clip is not None:
            try:
                expression_min_clip = float(expression_min_clip)
            except ValueError:
                return {
                    "success": -1,
                    "message": "The expression_min_clip parameter must be a number."
                }

        if not gene_symbol or not dataset_id:
            return {
                "success": -1,
                "message": "Request needs both dataset id and gene symbol."
            }

        dataset = geardb.get_dataset_by_id(dataset_id)

        if not dataset:
            return {
                "success": -1,
                "message": "The dataset was not found."
            }

        # ! SVG analysis does not operate on analysis objects.

        h5_path = dataset.get_file_path()
        if not os.path.exists(h5_path):
            return {
                "success": -1,
                "message": "The h5ad file was not found."
            }

        adata = sc.read_h5ad(h5_path)
        # convert adata.X to a dense matrix if it is sparse
        # This prevents issues with the min/max functions
        if scipy.sparse.issparse(adata.X):
            adata.X = adata.X.toarray() # type: ignore
        else:
            adata.X = np.asarray(adata.X)

        if projection_id:
            try:
                adata = create_projection_adata(adata, dataset_id, projection_id)
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }

        adata = clip_expression_values(adata, min_clip=expression_min_clip)

        gene_symbols = (gene_symbol,)

        if 'gene_symbol' not in adata.var.columns:
            return {"success": -1, "message": "The h5ad is missing the gene_symbol column."}

        gene_filter = adata.var.gene_symbol.isin(gene_symbols)
        if not gene_filter.any():
            return {
                "success": -1,
                "message": "The gene symbol '{}' was not found in the dataset.".format(gene_symbol)
            }

        # SAdkins - The code was rearranged because of how projection_adatas are handled
        # The projection_adata is file-backed to a tempfile and when the "selected" adata is created,
        # the tempfile is closed.  This causes the adata to be file-backed to a closed file.

        scores = {
            "dataset": dict()
            , "gene": dict()
            , "tissue": dict()
        }


        min_x = np.nanmin(adata.to_df().values)
        max_x = np.nanmax(adata.to_df().values)

        scores["dataset"] = {
            "dataset_id": dataset_id
            , "min": float(min_x)
            , "max": float(max_x)
        }

        tissues = adata.obs.index.tolist()
        for tissue in tissues:
            tissue_adata = adata[tissue, :]
            min_x = np.nanmin(tissue_adata.to_df().values)
            max_x = np.nanmax(tissue_adata.to_df().values)
            scores['tissue'][tissue] = {
                "min": float(min_x),
                "max": float(max_x)
            }

        try:
            selected = adata[:, gene_filter].to_memory()
        except Exception:
            # The "try" may fail for projections as it is already in memory
            print(f"Could not convert to memory for dataset {dataset_id}.  Using file-backed.", file=sys.stderr)
            selected = adata[:, gene_filter]

        df = selected.to_df()

        success = 1
        message = ""
        if len(df.columns) > 1:
            success = 2
            message = "WARNING: Multiple Ensemble IDs found for gene symbol '{}'.  Using the first stored Ensembl ID.".format(gene_symbol)
            df = df.iloc[:,[0]] # Note, put the '0' in a list to return a DataFrame.  Not having in list returns DataSeries instead

        scores["gene"] = {
                "gene": gene_symbol,
                "min": float(np.nanmin(selected.to_df().values)),
                "max": float(np.nanmax(selected.to_df().values))
            }

        if add_user_defined:
            # If either vmin or vmax is None, set to the min/max of the data
            # Also check if vmax is 0.0, which is the default value sent from the client
            if vmin is None or vmin == "0":
                vmin = 0.0
            if vmax is None or vmax == "0":
                vmax = scores["gene"]["max"]
            scores["user_defined"] = {
                "min": float(vmin),
                "max": float(vmax)
            }


        # Close adata so that we do not have a stale opened object
        if adata.isbacked:
            adata.file.close()

        # Get the average for all cells if there is a cell_type
        # in case there is an SVG path that has the convention
        # {tissue}--mean
        if 'cell_type' in selected.obs.columns:
            df = selected.to_df()
            df = pd.concat([df, selected.obs["cell_type"]], axis=1)

            cell_type_avgs = df.groupby('cell_type', observed=False).mean()
            # Add mean label
            mean_labels = cell_type_avgs.reset_index()['cell_type'].apply(
                lambda cell_type: f"{cell_type}--mean")
            cell_type_avgs = cell_type_avgs.set_index(mean_labels)

            # Remove column in order to match both df and cell_type_avgs
            # for concatenation
            del df['cell_type']
            # Concatenate our dataframe with cell_type averages with our
            # dataframe with all the tissues.
            df = pd.concat([df, cell_type_avgs])
        else:
            df = selected.to_df()
            df = pd.concat([df, selected.obs], axis=1)

        df = df.rename(columns={df.columns[0]: "data"})

        # Close adata so that we do not have a stale opened object
        if selected.isbacked:
            selected.file.close()

        return {
            "success": success,
            "message": message,
            "scores": scores,
            **df.to_dict()
        }
