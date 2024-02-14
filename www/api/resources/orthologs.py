from flask import request
from flask_restful import Resource
import os, sys
import geardb
from gear.orthology import get_ortholog_file, get_ortholog_files_from_dataset, map_single_gene, map_multiple_genes

def check_all_genes_in_dataset(adata, gene_symbols):
    """
    Check if all the given gene symbols are present in the dataset.

    Parameters:
    adata (AnnData): Annotated data object.
    gene_symbols (list): List of gene symbols to check.

    Returns:
    bool: True if all gene symbols are present, False otherwise.
    """
    gene_filter = adata.var.gene_symbol.isin(gene_symbols)
    return gene_filter.all()

def normalize_searched_gene(gene_list, chosen_gene):
    """Convert to case-insensitive version of gene.  Returns None if gene not found in dataset."""
    for g in gene_list:
        if chosen_gene.lower() == str(g).lower():
            return g
    return None

def get_mapped_gene_symbol(gene_symbol, gene_organism_id, dataset_organism_id):
    """
    Retrieves the mapped gene symbol for a given gene symbol, gene organism ID, and dataset organism ID.

    Args:
        gene_symbol (str): The gene symbol to be mapped.
        gene_organism_id (str): The organism ID of the gene.
        dataset_organism_id (str): The organism ID of the dataset.

    Returns:
        list: A list of mapped gene symbols.

    """
    if gene_organism_id and gene_organism_id != dataset_organism_id:
        ortholog_file = get_ortholog_file(gene_organism_id, dataset_organism_id, "ensembl")
        return map_single_gene(gene_symbol, ortholog_file)
    else:
        for ortholog_file in get_ortholog_files_from_dataset(dataset_organism_id, "ensembl"):
            try:
                mapped_genes = map_single_gene(gene_symbol, ortholog_file)
                if mapped_genes:
                    return mapped_genes
            except:
                continue
    return []

def get_mapped_gene_symbols(gene_symbols, gene_organism_id, dataset_organism_id):
    """
    Maps a list of gene symbols to their orthologous symbols in a given organism.

    Args:
        gene_symbols (list): List of gene symbols to be mapped.
        gene_organism_id (str): ID of the organism corresponding to the gene symbols.
        dataset_organism_id (str): ID of the organism corresponding to the dataset.

    Returns:
        dict: A dictionary mapping the input gene symbols to their orthologous symbols.
    """
    if gene_organism_id and gene_organism_id != dataset_organism_id:
        ortholog_file = get_ortholog_file(gene_organism_id, dataset_organism_id, "ensembl")

        return map_multiple_genes(gene_symbols, ortholog_file)
    else:
        for ortholog_file in get_ortholog_files_from_dataset(dataset_organism_id, "ensembl"):
            try:
                mapped_gene_symbols_dict =  map_multiple_genes(gene_symbols, ortholog_file)
                # ? Should we check all and return the dict with the most matches
                if len(mapped_gene_symbols_dict):
                    return mapped_gene_symbols_dict
            except Exception as e:
                print(str(e), file=sys.stderr)
                continue
    raise Exception("No orthologous mapping found for the given gene symbols. {}".format(gene_symbols))

def check_gene_in_dataset(adata, gene_symbol):
    """
    Check if any of the given gene symbols are present in the dataset.

    Args:
        adata (AnnData): Annotated data object.
        gene_symbols (list): List of gene symbols to check.

    Returns:
        bool: True if any of the gene symbols are present in the dataset, False otherwise.
    """

    dataset_genes = adata.var['gene_symbol'].unique().tolist()
    normalized_gene = normalize_searched_gene(dataset_genes, gene_symbol)

    gene_symbols = (normalized_gene,)
    gene_filter = adata.var.gene_symbol.isin(gene_symbols)
    return gene_filter.any()

class Orthologs(Resource):

    def get(self, dataset_id):
        gene_symbol = request.args.get('gene_symbol', None)
        gene_organism_id = request.args.get('gene_organism_id', None)

        if not gene_symbol:
            return {"error": "No gene symbol provided."}, 400

        if not dataset_id:
            return {"error": "No dataset ID provided."}, 400

        if gene_organism_id:
            gene_organism_id = int(gene_organism_id)

        # Get the dataset and organism ID
        dataset = geardb.get_dataset_by_id(dataset_id)
        dataset_organism_id = dataset.organism_id

        h5_path = dataset.get_file_path()
        if not os.path.exists(h5_path):
            return {"error": "The h5ad file was not found."}, 400

        import scanpy as sc
        adata = sc.read_h5ad(h5_path)

        dataset_genes = adata.var['gene_symbol'].unique().tolist()

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

        if check_gene_in_dataset(adata, gene_symbol):
            normalized_gene = normalize_gene(gene_symbol)
            mapped_gene_symbols_dict[gene_symbol] = [normalized_gene]
        else:
            try:
                mapped_gene_symbols = get_mapped_gene_symbol(gene_symbol, gene_organism_id, dataset_organism_id)
            except Exception as e:
                print(str(e), file=sys.stderr)
                return {"success": -1, "message": f"The searched gene symbol {gene_symbol} could not be mapped to the dataset organism."}

            # Filter out genes that are not in the dataset
            normalized_mapped_genes = [normalize_gene(mapped_gene_symbol) for mapped_gene_symbol in mapped_gene_symbols if check_gene_in_dataset(adata, mapped_gene_symbol)]
            mapped_gene_symbols_dict[gene_symbol] = normalized_mapped_genes

        return {"success": 1, "mapping": mapped_gene_symbols_dict}, 200


    def post(self, dataset_id):
        session_id = request.cookies.get('gear_session_id')
        req = request.get_json()
        gene_symbols = req.get('gene_symbols', None)
        analysis = req.get('analysis', None)
        gene_organism_id = req.get('gene_organism_id', None)

        if not gene_symbols:
            return {"error": "No gene symbols provided."}, 400

        if not dataset_id:
            return {"error": "No dataset ID provided."}, 400

        if gene_organism_id:
            gene_organism_id = int(gene_organism_id)

        # Get the dataset and organism ID
        dataset = geardb.get_dataset_by_id(dataset_id)
        dataset_organism_id = dataset.organism_id

        # Get the right AnnData object depending on if the analysis is provided
        try:
            ana = geardb.get_analysis(analysis, dataset_id, session_id)
        except Exception as e:
            return {"error": str(e)}, 400

        adata = ana.get_adata(backed=False)
        dataset_genes = adata.var['gene_symbol'].unique().tolist()

        def normalize_gene(gene):
            return normalize_searched_gene(dataset_genes, gene)

        mapped_gene_symbols_dict = {gene_symbol: [] for gene_symbol in gene_symbols}

        if gene_organism_id and gene_organism_id == dataset_organism_id:
            # Using adata with "backed" mode does not work with volcano plot
            for gene_symbol in gene_symbols:
                normalized_gene = normalize_gene(gene_symbol)
                if normalized_gene:
                    mapped_gene_symbols_dict[gene_symbol] = normalized_gene
        else:
            # Perform orthology mapping
            mapped_gene_symbols_dict = get_mapped_gene_symbols(gene_symbols, gene_organism_id, dataset_organism_id)

            # for each mapped gene symbol, verify the mapped genes are in the dataset and normalize to those genes
            for gene_symbol in gene_symbols:
                normalized_mapped_genes = [normalize_gene(mapped_gene_symbol) for mapped_gene_symbol in mapped_gene_symbols_dict[gene_symbol] if check_gene_in_dataset(adata, mapped_gene_symbol)]
                mapped_gene_symbols_dict[gene_symbol] = normalized_mapped_genes

        # Return a dictionary where the key is the original gene symbol name and the mapping is a list of orthologs
        return {"success": 1, "mapping": mapped_gene_symbols_dict}, 200
