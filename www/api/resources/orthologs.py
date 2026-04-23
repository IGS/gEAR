import sys

import geardb
from flask import request
from flask_restful import Resource
from gear.orthology import (
    get_ortholog_file,
    get_ortholog_files_from_dataset,
    map_multiple_genes,
    map_single_gene,
)
from gear.utils import catch_memory_error

from .common import get_adata_shadow, get_spatial_adata

def fetch_ortholog_files(gene_organism_id, dataset_organism_id, exclusive_org=False):
    """
    Fetch ortholog files based on organism IDs.

    Prioritizes the gene organism's ortholog file if gene_organism_id is provided.

    Args:
        gene_organism_id (int): Organism ID of the gene (optional).
        dataset_organism_id (int): Organism ID of the dataset.
        exclusive_org (bool): If True, only fetch from gene organism.

    Returns:
        list: List of ortholog file paths, with gene organism file first if available.

    Raises:
        FileNotFoundError: If ortholog files cannot be found.
    """
    is_cross_organism = gene_organism_id and gene_organism_id != dataset_organism_id
    ortholog_files = []
    try:
        # Prioritize gene organism file if cross-organism
        if is_cross_organism:
            gene_org_file = get_ortholog_file(gene_organism_id, dataset_organism_id, "ensembl")
            if gene_org_file:
                ortholog_files.append(gene_org_file)

        # Add dataset organism files unless exclusive
        if not exclusive_org:
            ortholog_files.extend(get_ortholog_files_from_dataset(dataset_organism_id, "ensembl"))
        return ortholog_files
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        raise

def get_mapped_gene_symbols(gene_symbols, gene_organism_id, dataset_organism_id, exclusive_org=False):
    """Retrieves the mapped gene symbols for the given gene symbols."""
    try:
        ortholog_files = fetch_ortholog_files(gene_organism_id, dataset_organism_id, exclusive_org)
    except FileNotFoundError:
        return {gene: [] for gene in gene_symbols}

    for ortholog_file in ortholog_files:
        mapped_genes_dict = map_multiple_genes(gene_symbols, ortholog_file)

        if any(len(v) > 0 for v in mapped_genes_dict.values()):
            return mapped_genes_dict

        if gene_organism_id and exclusive_org:
            raise Exception(f"No orthologous mapping found for gene symbols: {gene_symbols}.")

    return {gene: [] for gene in gene_symbols}

def build_gene_map(dataset_genes: set) -> dict:
    """
    Build normalized gene map with helper function.

    Returns:
        dict: Contains 'gene_map_set' and 'normalize_gene' callable.
    """
    gene_map = {str(g).lower(): g for g in dataset_genes}
    return {
        'gene_map_set': set(gene_map.keys()),
        'normalize_gene': lambda gene: gene_map.get(str(gene).lower())
    }

def normalize_mapped_genes(mapped_gene_symbols: list, gene_set: set, normalize_gene) -> list:
    """
    Filter and normalize mapped gene symbols to those present in the dataset.

    Args:
        mapped_gene_symbols (list): List of gene symbols that were mapped via orthology.
        gene_set (set): Set of lowercase gene symbols present in the dataset.
        normalize_gene (callable): Function to normalize a gene symbol using gene_set.

    Returns:
        list: Normalized gene symbols that exist in the dataset.
    """
    normalized_genes = []
    for mapped_gene_symbol in mapped_gene_symbols:
        if check_gene_in_dataset(gene_set, mapped_gene_symbol):
            normalized_gene = normalize_gene(mapped_gene_symbol)
            if normalized_gene is not None:
                normalized_genes.append(normalized_gene)
    return normalized_genes

def check_gene_in_dataset(gene_map: set, gene_symbol: str) -> bool:
    """
    Check whether a gene symbol is present in a dataset gene map.

    Parameters
    ----------
    gene_map : set
        A set containing gene symbols (expected to be lowercase strings) that
        represent the genes present in the dataset. Membership is tested
        directly against this set.
    gene_symbol : str
        The gene symbol to check. This value will be coerced to str and
        lowercased before testing membership.

    Returns
    -------
    bool
        True if the lowercased gene_symbol is found in gene_map, False otherwise.

    Notes
    -----
    - Non-string inputs for gene_symbol are converted via str() before lowercasing.
    - For correct case-insensitive behavior, gene_map should contain lowercase
      representations of gene symbols.
    - This function performs a simple membership test and does not perform
      additional normalization (e.g., trimming whitespace, handling synonyms).

    Examples
    --------
    >>> check_gene_in_dataset({'tp53', 'brca1'}, 'TP53')
    True
    >>> check_gene_in_dataset({'tp53', 'brca1'}, 123)
    False
    """
    return str(gene_symbol).lower() in gene_map


def load_dataset_genes(analysis_id, dataset_id, session_id, dtype):
    """Load genes from dataset file."""
    import gc

    if dtype == "spatial":
        adata = get_spatial_adata(analysis_id, dataset_id, session_id)
    else:
        adata = get_adata_shadow(analysis_id, dataset_id, session_id)

    if not hasattr(adata, 'var') or 'gene_symbol' not in adata.var:
        raise ValueError("Dataset does not contain 'gene_symbol' field.")

    dataset_genes = set(adata.var['gene_symbol'].unique())

    try:
        adata.file.close()
    except AttributeError:
        del adata
        gc.collect()

    return dataset_genes


def map_all_genes(gene_symbols: list, gene_organism_id: int, dataset_organism_id: int, gene_map_set: set, exclusive_org: bool, normalize_gene):
    """Map all genes using direct and orthology methods."""
    mapped_dict = {symbol: [] for symbol in gene_symbols}
    unmapped = []

    # Direct mapping (from dataset)
    if gene_organism_id and gene_organism_id == dataset_organism_id:
        for symbol in gene_symbols:
            normalized = normalize_gene(symbol)
            if normalized:
                mapped_dict[symbol] = [normalized]
            else:
                unmapped.append(symbol)

        if not unmapped:
            return mapped_dict
    else:
        unmapped = gene_symbols

    # Orthology mapping
    try:
        ortho_mapped = get_mapped_gene_symbols(unmapped, gene_organism_id, dataset_organism_id, exclusive_org)
        for symbol, genes in ortho_mapped.items():
            # Attempt to normalize mapped genes and filter to those in the dataset. If none remain, try last chance mapping.
            normalized = normalize_mapped_genes(genes, gene_map_set, normalize_gene)
            mapped_dict[symbol] = normalized if normalized else last_chance_map(symbol, gene_map_set, normalize_gene)
    except Exception as e:
        print(str(e), file=sys.stderr)
        if exclusive_org:
            raise

    return mapped_dict

def last_chance_map(gene_symbol: str, gene_map_set: set, normalize_gene):
    """Final attempt to map gene if orthology failed."""
    if check_gene_in_dataset(gene_map_set, gene_symbol):
        return [normalize_gene(gene_symbol)]
    return []

class Orthologs(Resource):

    @catch_memory_error()
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
        if not req:
            return {"error": "Invalid JSON body."}, 400
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

        if not dataset.organism_id:
            return {"error": "The dataset does not have an associated organism."}, 400

        try:
            dataset_genes = load_dataset_genes(analysis_id, dataset_id, session_id, dataset.dtype)
        except FileNotFoundError:
            return {"success": -1, "message": "No dataset file found."}, 400
        except Exception as e:
            return {"success": -1, "message": str(e)}, 400

        # Map genes
        gene_map_data = build_gene_map(dataset_genes)
        result = map_all_genes(
            gene_symbols,
            gene_organism_id,
            dataset.organism_id,
            gene_map_data['gene_map_set'],
            exclusive_org,
            gene_map_data['normalize_gene'],
        )

        return {"success": 1, "mapping": result}, 200
