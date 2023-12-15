import os
from pathlib import Path

import geardb
import numpy as np
import pandas as pd
import scanpy as sc
from flask import request
from flask_restful import Resource

from werkzeug.utils import secure_filename

from gear.orthology import get_ortholog_file, get_ortholog_files_from_dataset, map_single_gene

TWO_LEVELS_UP = 2
abs_path_www = Path(__file__).resolve().parents[TWO_LEVELS_UP] # web-root dir
PROJECTIONS_BASE_DIR = abs_path_www.joinpath('projections')

class PlotError(Exception):
    """Error based on plotting issues."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)

def get_mapped_gene_symbol(gene_symbol, gene_organism_id, dataset_organism_id):
    """
    Maps a gene symbol to its corresponding orthologous gene symbol in a given dataset.

    Args:
        gene_symbol (str): The gene symbol to be mapped.
        gene_organism_id (str): The organism ID of the gene symbol.
        dataset_organism_id (str): The organism ID of the dataset.

    Returns:
        str: The mapped orthologous gene symbol, or None if no mapping is found.
    """
    if gene_organism_id and gene_organism_id != dataset_organism_id:
        ortholog_file = get_ortholog_file(gene_organism_id, dataset_organism_id, "ensembl")
        return map_single_gene(gene_symbol, ortholog_file)
    else:
        for ortholog_file in get_ortholog_files_from_dataset(dataset_organism_id, "ensembl"):
            try:
                return map_single_gene(gene_symbol, ortholog_file)
            except:
                continue
    return None

def check_gene_in_dataset(adata, gene_symbols):
    """
    Check if any of the given gene symbols are present in the dataset.

    Args:
        adata (AnnData): Annotated data object.
        gene_symbols (list): List of gene symbols to check.

    Returns:
        bool: True if any of the gene symbols are present in the dataset, False otherwise.
    """
    gene_filter = adata.var.gene_symbol.isin(gene_symbols)
    return gene_filter.any()

def create_projection_adata(dataset_adata, dataset_id, projection_id):
    # Create AnnData object out of readable CSV file
    # ? Does it make sense to put this in the geardb/Analysis class?
    projection_id = secure_filename(projection_id)
    dataset_id = secure_filename(dataset_id)

    projection_dir = Path(PROJECTIONS_BASE_DIR).joinpath("by_dataset", dataset_id)
    # Sanitize input to prevent path traversal
    projection_adata_path = projection_dir.joinpath("{}.h5ad".format(projection_id))

    if projection_adata_path.is_file():
        return sc.read_h5ad(projection_adata_path)  #, backed="r")

    projection_csv_path = projection_dir.joinpath("{}.csv".format(projection_id))
    try:
        projection_adata = sc.read_csv(projection_csv_path)
    except:
        raise PlotError("Could not create projection AnnData object from CSV.")
    projection_adata.obs = dataset_adata.obs
    # Close dataset adata so that we do not have a stale opened object
    if dataset_adata.isbacked:
        dataset_adata.file.close()
    projection_adata.var["gene_symbol"] = projection_adata.var_names
    # Associate with a filename to ensure AnnData is read in "backed" mode
    projection_adata.filename = projection_adata_path
    return projection_adata

class SvgData(Resource):
    """Resource for retrieving data from h5ad to be used to color svgs.

    Returns
    -------
    dict
        SVG data
    """
    def get(self, dataset_id):
        gene_symbol = request.args.get('gene', None)
        gene_organism_id = request.get('gene_organism_id', None)
        projection_id = request.args.get('projection_id', None)    # projection id of csv output

        if not gene_symbol or not dataset_id:
            return {
                "success": -1,
                "message": "Request needs both dataset id and gene symbol."
            }

        dataset = geardb.get_dataset_by_id(dataset_id)
        dataset_organism_id = dataset.organism_id

        h5_path = dataset.get_file_path()
        if not os.path.exists(h5_path):
            return {
                "success": -1,
                "message": "The h5ad file was not found."
            }

        adata = sc.read_h5ad(h5_path)

        if projection_id:
            try:
                adata = create_projection_adata(adata, dataset_id, projection_id)
            except PlotError as pe:
                return {
                    'success': -1,
                    'message': str(pe),
                }


        mapped_gene_symbol = None
        gene_symbols = (gene_symbol,)

        if 'gene_symbol' not in adata.var.columns:
            return {"success": -1, "message": "The h5ad is missing the gene_symbol column."}

        if not check_gene_in_dataset(adata, gene_symbols):
            try:
                mapped_gene_symbol = get_mapped_gene_symbol(gene_symbol, gene_organism_id, dataset_organism_id)
            except:
                return {"success": -1, "message": f"The searched gene symbol {gene_symbol} could not be mapped to the dataset organism."}

            if mapped_gene_symbol:
                gene_symbols = (mapped_gene_symbol,)
                if not check_gene_in_dataset(adata, gene_symbols):
                    return {"success": -1, "message": f"The searched gene symbol {gene_symbol} could not be found in the h5ad file."}
            else:
                return {"success": -1, "message": f"The searched gene symbol {gene_symbol} could not be mapped to the dataset organism."}

        try:
            gene_filter = adata.var.gene_symbol.isin(gene_symbols)
            selected = adata[:, gene_filter].to_memory()
        except:
            # The "try" may fail for projections as it is already in memory
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

        # Close adata so that we do not have a stale opened object
        if adata.isbacked:
            adata.file.close()

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

        # Close adata so that we do not have a stale opened object
        if selected.isbacked:
            selected.file.close()

        return {
            "success": success,
            "message": message,
            "scores": scores,
            "mapped_gene_symbol": mapped_gene_symbol,
            **df.to_dict()
        }
