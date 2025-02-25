from flask import request
from flask_restful import Resource
import os, sys
import geardb
from gear.orthology import get_ortholog_file, get_ortholog_files_from_dataset, map_single_gene, map_multiple_genes
from .common import get_adata_shadow

def normalize_searched_gene(gene_set, chosen_gene):
    """Convert to case-insensitive version of gene.  Returns None if gene not found in dataset."""
    # Some genes get read in as floats by accident
    chosen_gene_lower = str(chosen_gene).lower()
    for gene in gene_set:
        try:
            if chosen_gene_lower == str(gene).lower():
                return gene
        except Exception:
            print(gene, file=sys.stderr)
            raise
    return None

def get_mapped_gene_symbol(gene_symbol, gene_organism_id, dataset_organism_id, exclusive_org=False):
    """
    Retrieves the mapped gene symbol for a given gene symbol, gene organism ID, and dataset organism ID.

    Args:
        gene_symbol (str): The gene symbol to be mapped.
        gene_organism_id (str): The organism ID of the gene.
        dataset_organism_id (str): The organism ID of the dataset.
        exclusive_org (bool, optional): Flag indicating whether to only consider orthologs from the gene organism.
            Defaults to False.

    Returns:
        list or dict: The mapped gene symbol(s) if found, otherwise an empty list or dict.

    Raises:
        Exception: If no orthologous mapping is found for the given gene symbol.
    """

    def fetch_ortholog_files():
        # Determine if we need to get a single ortholog file or multiple
        if gene_organism_id and gene_organism_id != dataset_organism_id:
            # Get a single ortholog file
            ortholog_files = [get_ortholog_file(gene_organism_id, dataset_organism_id, "ensembl")]
            if not exclusive_org:
                ortholog_files += get_ortholog_files_from_dataset(dataset_organism_id, "ensembl")
        else:
            ortholog_files = get_ortholog_files_from_dataset(dataset_organism_id, "ensembl")
        return ortholog_files

    try:
        ortholog_files = fetch_ortholog_files()
    except FileNotFoundError as e:
        # We want this to fail gracefully, so return an empty list. The original will be mapped back to itself downstream.
        return []

    for ortholog_file in ortholog_files:
        mapped_genes = map_single_gene(gene_symbol, ortholog_file)
        if mapped_genes:
            return mapped_genes

        # At this point, we have an empty list or dict, so we should continue to the next ortholog file
        if gene_organism_id and exclusive_org:
            raise Exception(f"No orthologous mapping found for the given gene symbols {gene_symbol}.")

    return []


def get_mapped_gene_symbols(gene_symbols, gene_organism_id, dataset_organism_id, exclusive_org=False):
    """
    Retrieves the mapped gene symbols for the given gene symbols, gene organism ID, and dataset organism ID.

    Args:
        gene_symbols (list): List of gene symbols to map.
        gene_organism_id (str): Organism ID of the gene.
        dataset_organism_id (str): Organism ID of the dataset.
        exclusive_org (bool, optional): Flag indicating whether to only consider orthologs from the gene organism.
            Defaults to False.

    Returns:
        dict: A dictionary containing the mapped gene symbols as keys and their corresponding orthologs as values.

    Raises:
        Exception: If no orthologous mapping is found for the given gene symbols.
    """
    # Determine if we need to get a single ortholog file or multiple
    is_single_ortholog_file_needed = gene_organism_id and gene_organism_id != dataset_organism_id

    try:
        if is_single_ortholog_file_needed:
            # Get a single ortholog file
            ortholog_files = [get_ortholog_file(gene_organism_id, dataset_organism_id, "ensembl")]
            if not exclusive_org:
                ortholog_files += get_ortholog_files_from_dataset(dataset_organism_id, "ensembl")

        else:
            # Get multiple ortholog files from the dataset
            ortholog_files = get_ortholog_files_from_dataset(dataset_organism_id, "ensembl")
    except FileNotFoundError as e:
        # We want this to fail gracefully, so return an empty list for each gene symbol
        # The original will be mapped back to itself downstream.
        return {gene: [] for gene in gene_symbols}

    for ortholog_file in ortholog_files:
        mapped_genes_dict = map_multiple_genes(gene_symbols, ortholog_file)

        # If any keys have a non-empty list, return the dict
        if any([len(v) > 0 for v in mapped_genes_dict.values()]):
            return mapped_genes_dict
        # At this point, we have an empty list or dict, so we should continue to the next ortholog file
        if gene_organism_id and exclusive_org:
            raise Exception(f"No orthologous mapping found for the given gene symbols {gene_symbols}.")
        continue

    return {gene: [] for gene in gene_symbols}


def check_gene_in_dataset(adata, gene_symbol):
    """
    Check if any of the given gene symbols are present in the dataset.

    Args:
        adata (AnnData): Annotated data object.
        gene_symbols (list): List of gene symbols to check.

    Returns:
        bool: True if any of the gene symbols are present in the dataset, False otherwise.
    """

    dataset_genes = set(adata.var['gene_symbol'].unique())
    normalized_gene = normalize_searched_gene(dataset_genes, gene_symbol)

    gene_symbols = (normalized_gene,)
    gene_filter = adata.var.gene_symbol.isin(gene_symbols)
    return gene_filter.any()

class Orthologs(Resource):

    def get(self, dataset_id):
        """
        Retrieve gene symbol mappings for a given dataset. This is meant to search a single gene symbol

        Args:
            dataset_id (str): The ID of the dataset to query.

        Returns:
            tuple: A tuple containing a dictionary with the result and an HTTP status code.
            - If successful, returns a dictionary with the key "success" set to 1 and the key "mapping" containing the mapped gene symbols.
            - If an error occurs, returns a dictionary with the key "error" or "message" and an appropriate HTTP status code.

        Raises:
            Exception: If an error occurs during orthology mapping and exclusive_org is True.

        Notes:
            - The function expects the following query parameters:
            - gene_symbol (str): The gene symbol to search for.
            - gene_organism_id (int, optional): The organism ID of the gene symbol.
            - exclusive_org (str, optional): If "true", return an error if mapping is not found. Defaults to "false".
        """
        gene_symbol = request.args.get('gene_symbol')
        gene_organism_id = request.args.get('gene_organism_id' )
        exclusive_org = request.args.get('exclusive_org', "false").lower() == "true" # If true, only check this organism. If false, loop through other organisms

        if not gene_symbol:
            return {"error": "No gene symbol provided."}, 400

        if not dataset_id:
            return {"error": "No dataset ID provided."}, 400

        if gene_organism_id:
            gene_organism_id = int(gene_organism_id)


        # Get the dataset and organism ID
        dataset = geardb.get_dataset_by_id(dataset_id)
        if not dataset:
            return {"error": "The dataset was not found."}, 400
        dataset_organism_id = dataset.organism_id

        h5_path = dataset.get_file_path()
        if not os.path.exists(h5_path):
            return {"error": "The h5ad file was not found."}, 400

        import scanpy as sc
        adata = sc.read_h5ad(h5_path)

        dataset_genes = set(adata.var['gene_symbol'].unique())

        def normalize_gene(gene):
            return normalize_searched_gene(dataset_genes, gene)

        mapped_gene_symbols_dict = {gene_symbol: []}

        if gene_organism_id and gene_organism_id == dataset_organism_id:
            normalized_gene = normalize_gene(gene_symbol)
            if normalized_gene:
                mapped_gene_symbols_dict[gene_symbol] = [normalized_gene]
                return {"success": 1, "mapping": mapped_gene_symbols_dict}, 200
            else:
                return {"success": -1, "message": f"The searched gene symbol {gene_symbol} could not be found in the dataset."}

        # Perform orthology mapping.
        try:
            mapped_gene_symbols = get_mapped_gene_symbol(gene_symbol, gene_organism_id, dataset_organism_id, exclusive_org)
        except Exception as e:
            if exclusive_org:
                return {"success": -1, "message": str(e)}
            else:
                return {"error": str(e)}, 400

        # Filter out genes that are not in the dataset
        normalized_mapped_genes = [normalize_gene(mapped_gene_symbol) for mapped_gene_symbol in mapped_gene_symbols if check_gene_in_dataset(adata, mapped_gene_symbol)]
        mapped_gene_symbols_dict[gene_symbol] = normalized_mapped_genes

        # last chance to map.  Check if nonmapping genes are actually in the dataset (since gene_organism_id may not have been provided)
        if not normalized_mapped_genes:
            if check_gene_in_dataset(adata, gene_symbol):
                normalized_gene = normalize_gene(gene_symbol)
                mapped_gene_symbols_dict[gene_symbol] = [normalized_gene]
            else:
                return {"success": -1, "message": f"The searched gene symbol {gene_symbol} could not be found in the dataset."}

        return {"success": 1, "mapping": mapped_gene_symbols_dict}, 200


    def post(self, dataset_id):
        """
        Handles POST requests to map gene symbols to their orthologs within a specified dataset.

        Args:
            dataset_id (str): The ID of the dataset to be queried.

        Returns:
            tuple: A tuple containing a JSON response and an HTTP status code.
                - On success, returns a JSON object with a "success" key set to 1 and a "mapping" key containing a dictionary
                  where the keys are the original gene symbols and the values are lists of mapped orthologs.
                - On failure, returns a JSON object with an "error" key and an appropriate error message, along with a 400 status code.

        Request JSON Structure:
            {
                "gene_symbols": list of str,  # List of gene symbols to be mapped.
                "analysis": str,              # Analysis identifier (optional).
                "gene_organism_id": int,      # Organism ID for the genes (optional).
                "exclusive_org": str          # "true" or "false" indicating whether to exclusively check the specified organism (default: "false").
            }

        Raises:
            Exception: If there is an error in retrieving the analysis or performing the orthology mapping.
        """
        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        gene_symbols = req.get('gene_symbols')
        analysis = req.get('analysis')
        gene_organism_id = req.get('gene_organism_id')
        exclusive_org = req.get('exclusive_org', "false").lower() == "true" # If true, only check this organism. If false, loop through other organisms

        if not gene_symbols:
            return {"error": "No gene symbols provided."}, 400

        if not dataset_id:
            return {"error": "No dataset ID provided."}, 400

        if gene_organism_id:
            gene_organism_id = int(gene_organism_id)

        analysis_id = None
        if analysis:
            analysis_id = analysis.get('id')

        # Get the dataset and organism ID
        dataset = geardb.get_dataset_by_id(dataset_id)

        if not dataset:
            return {"error": "The dataset was not found."}, 400

        dataset_path = dataset.get_file_path()
        if not os.path.exists(dataset_path):
            return {"error": f"The dataset for id {dataset_id} was not found."}, 400

        dataset_organism_id = dataset.organism_id

        try:
            adata = get_adata_shadow(analysis_id, dataset_id, session_id, dataset_path, include_images=False)
        except FileNotFoundError:
            return {
                "success": -1,
                'message': "No h5 file found for this dataset"
            }

        dataset_genes = set(adata.var['gene_symbol'].unique())

        def normalize_gene(gene):
            return normalize_searched_gene(dataset_genes, gene)

        mapped_gene_symbols_dict = {gene_symbol: [] for gene_symbol in gene_symbols}

        if gene_organism_id and gene_organism_id == dataset_organism_id:
            # Using adata with "backed" mode does not work with volcano plot
            for gene_symbol in gene_symbols:
                normalized_gene = normalize_gene(gene_symbol)
                if normalized_gene:
                    mapped_gene_symbols_dict[gene_symbol] = [normalized_gene]

            return {"success": 1, "mapping": mapped_gene_symbols_dict}, 200

        # Perform orthology mapping.
        try:
            mapped_gene_symbols_dict = get_mapped_gene_symbols(gene_symbols, gene_organism_id, dataset_organism_id, exclusive_org)
        except Exception as e:
            if exclusive_org:
                return {"success": -1, "message": str(e)}
            else:
                return {"error": str(e)}, 400

        genes_not_mapped = []

        # for each mapped gene symbol, verify the mapped genes are in the dataset and normalize to those genes
        for gene_symbol in gene_symbols:
            normalized_mapped_genes = [normalize_gene(mapped_gene_symbol) for mapped_gene_symbol in mapped_gene_symbols_dict[gene_symbol] if check_gene_in_dataset(adata, mapped_gene_symbol)]

            if not normalized_mapped_genes:
                genes_not_mapped.append(gene_symbol)
                continue

            mapped_gene_symbols_dict[gene_symbol] = normalized_mapped_genes

        # last chance to map.  Check if nonmapping genes are actually in the dataset (since gene_organism_id may not have been provided)
        for gene_symbol in genes_not_mapped:
            if check_gene_in_dataset(adata, gene_symbol):
                normalized_gene = normalize_gene(gene_symbol)
                mapped_gene_symbols_dict[gene_symbol] = [normalized_gene]
            else:
                mapped_gene_symbols_dict[gene_symbol] = []

        # Return a dictionary where the key is the original gene symbol name and the mapping is a list of orthologs
        return {"success": 1, "mapping": mapped_gene_symbols_dict}, 200
