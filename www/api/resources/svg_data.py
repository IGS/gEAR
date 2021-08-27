from flask import request, make_response, jsonify, abort
from flask_restful import Resource
import numpy as np
import pandas as pd
import scanpy as sc
import json
import os
import sys
import geardb

class SvgData(Resource):
    """Resource for retrieving data from h5ad to be used to color svgs.

    Returns
    -------
    dict
        SVG data
    """
    def get(self, dataset_id):
        gene_symbol = request.args.get('gene')
        if not gene_symbol or not dataset_id:
            return {
                "success": -1,
                "message": "Request needs both dataset id and gene symbol."
            }

        dataset = geardb.get_dataset_by_id(dataset_id)
        # If the dataset is not public, make sure the user
        # requesting resource owns the dataset
        # This is commented right now until we work out modeling URL-shared datasets/profiles
        #if not dataset.is_public:
        #    session_id = request.cookies.get('gear_session_id')
        #    user = geardb.get_user_from_session_id(session_id)
        #    if user.id != dataset.owner_id:
        #        return {
        #            "success": -1,
        #            "message": 'Only the owner can access this dataset.'
        #        }


        h5_path = dataset.get_file_path()
        if not os.path.exists(h5_path):
            return {
                "success": -1,
                "message": "The h5ad file was not found."
            }

        adata = sc.read_h5ad(h5_path)
        gene_symbols = (gene_symbol,)

        if 'gene_symbol' in adata.var.columns:
            gene_filter = adata.var.gene_symbol.isin(gene_symbols)
            if not gene_filter.any():
                return {
                    "success": -1,
                    "message": "Gene not found."
                }
        else:
            return {
                "success": -1,
                "message": "The h5ad is missing gene_symbol.",
            }

        selected = adata[:, gene_filter]

        df = selected.to_df()

        success = 1
        message = ""
        if len(df.columns) > 1:
            success = 2
            message = "WARNING: Multiple Ensemble IDs found for gene symbol '{}'.  Using the first stored Ensembl ID.".format(gene_symbol)
            df = df.iloc[:,[0]] # Note, put the '0' in a list to return a DataFrame.  Not having in list returns DataSeries instead

        scores = {
          "dataset": {
            "dataset_id": dataset_id,
            "min": float(adata.X[~np.isnan(adata.X)].min()),
            "max": float(adata.X[~np.isnan(adata.X)].max())
          },
          "gene": {
            "gene": gene_symbol,
            "min": float(selected.X[~np.isnan(selected.X)].min()),
            "max": float(selected.X[~np.isnan(selected.X)].max())
          },
          "tissue": dict()
        }

        tissues = adata.obs.index.tolist()
        for tissue in tissues:
          tissue_adata = adata[tissue, :]
          scores['tissue'][tissue] = {
            "min": float(tissue_adata.X[~np.isnan(tissue_adata.X)].min()),
            "max": float(tissue_adata.X[~np.isnan(tissue_adata.X)].max())
          }

        # Get the average for all cells if there is a cell_type
        # in case there is an SVG path that has the convention
        # {tissue}--mean
        if 'cell_type' in selected.obs.columns:
          df = selected.to_df()
          df = pd.concat([df, selected.obs["cell_type"]], axis=1)

          cell_type_avgs = df.groupby('cell_type').mean()
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

        return {
          "success": success,
          "message": message,
          "scores": scores,
          **df.to_dict()
        }
