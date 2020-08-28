from flask import request
from flask_restful import Resource
import pandas as pd
import scanpy as sc
import json
import os
import sys
import copy
import geardb


class EpivizData(Resource):
    """Resource for retrieving data from h5ad to be used to draw charts on UI.
    Parameters
    ----------
    gene_symbol: str
        Gene to search in adata.
    index: list[str]
        List of column names to index on.
    colors: dict
        Dictionary mapping levels to colors.
    order: dict
        Dictionary mapping order levels.
    plot_type: str
        Plot type (bar, violin, scatter or line)
    Returns
    -------
    dict
        Plot data
    """

    def get(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        gene_symbol = request.args.get('gene')
        genome = request.args.get('genome')

        # TODO: may be loading all genomes before hand is faster ?
        # also needs a fallback is the gene name cannot be found
        genes =  pd.read_csv("/var/www/epiviz-api/genomes/" + genome + "/" + genome + ".txt", sep="\t",
            names=["chr", "start", "end", "width", "strand", "geneid", "exon_starts", "exon_ends", "gene"])
        
        if len(gene_symbol) > 1:
            genes = genes[genes['gene'].str.contains(gene_symbol, na=False, case=False)]

            if len(genes) > 0:
                return {
                    "chr": genes.loc[genes.index[0], 'chr'],
                    "start": int(genes.loc[genes.index[0], 'start']),
                    "end": int(genes.loc[genes.index[0], 'end'])
                }
            else :
                return {
                    "success": -1,
                    "message": "cannot find any gene symbol match in " + genomes
                }
        else :
            return {
                    "success": -1,
                    "message": "cannot find any gene symbol match in " + genomes
                }
