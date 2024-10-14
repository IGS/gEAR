# orthology.py - functions related to orthology mapping

"""
This is a library of functions related to orthology mapping. It is used to map gene symbols to their orthologous gene symbols in different organisms.
The ortholog mappings are stored in HDF5 files, and contain the following columns:
- id1: Ensembl ID of the gene in the first organism
- gs1: Gene symbol of the gene in the first organism
- id2: Ensembl ID of the gene in the second organism
- gs2: Gene symbol of the gene in the second organism
- algorithms_match_count: The number of algorithms that matched the orthologous gene symbols. Used to rank the orthologous gene symbols.
"""


import sys
import pandas as pd

from pathlib import Path

# append parent directory to path
sys.path.append(str(Path(__file__).resolve().parents[1]))

import geardb

organisms = geardb.OrganismCollection().get_all()

# Calculate these paths once at the module level
abs_path_www = Path(__file__).resolve().parents[2].joinpath("www") # get www directory
ORTHOLOG_BASE_DIR = abs_path_www.joinpath("feature_mapping")

def format_orthomap_file_base(first_org_id, second_org_id, annotation_source):
    """
    Formats the base filename for an orthomap file.

    Args:
        first_org_id (str): The ID of the first organism.
        second_org_id (str): The ID of the second organism.
        annotation_source (str): The source of the annotation.

    Returns:
        str: The formatted base filename for the orthomap file.
    """
    return "orthomap.{0}.{2}__{1}.{2}.hdf5".format(first_org_id, second_org_id, annotation_source)

def get_ortholog_file(gene_organism_id: str, dataset_organism_id: str, annotation_source: str="ensembl"):
    """
    Get the ortholog file for a given gene organism ID, dataset organism ID, and annotation source.

    Args:
        gene_organism_id (str): The ID of the gene organism.
        dataset_organism_id (str): The ID of the dataset organism.
        annotation_source (str, optional): The annotation source. Defaults to "ensembl".

    Returns:
        pathlib.Path: The path to the ortholog file.

    Raises:
        FileNotFoundError: If the orthologous mapping file is not found.
    """
    orthomap_file_base = format_orthomap_file_base(gene_organism_id, dataset_organism_id, annotation_source)
    orthomap_file = ORTHOLOG_BASE_DIR / orthomap_file_base

    if not orthomap_file.is_file():
        raise FileNotFoundError(f"Orthologous mapping file not found: {orthomap_file}")
    return orthomap_file

def get_ortholog_files_from_dataset(dataset_organism_id: str, annotation_source: str="ensembl"):
    """
    Get a list of orthologous mapping files for a given dataset organism ID and annotation source.

    Args:
        dataset_organism_id (str): The ID of the dataset organism.
        annotation_source (str, optional): The annotation source. Defaults to "ensembl".

    Returns:
        list: A list of tuples, where each tuple contains the path to an orthologous mapping file and its inverse.

    Raises:
        FileNotFoundError: If no orthologous mapping files are found for the dataset organism.
    """
    orthomap_file_base = format_orthomap_file_base( "[0-9]", dataset_organism_id, annotation_source)
    orthomap_files = list(ORTHOLOG_BASE_DIR.glob(orthomap_file_base))

    # sort by gene organism id
    orthomap_files.sort(key=lambda x: int(x.name.split("__")[0].split(".")[1]))

    if not orthomap_files:
        dataset_organism_name = get_organism_name_by_id(dataset_organism_id)
        raise FileNotFoundError(f"Orthologous mapping files not found for dataset organism {dataset_organism_name}")

    return list(orthomap_files)

def filter_organism_by_id(organism_id: str):
    """Filter the organisms list by the given organism ID.

    Args:
        organism_id (str): The organism ID to filter by.

    Returns:
        dict: The organism dictionary corresponding to the given organism ID.
    """

    return next((item for item in organisms if item.id == organism_id), None)

def get_organism_name_by_id(organism_id: str):
    """Get the organism name corresponding to the given organism ID.

    Args:
        organism_id (str): The organism ID to filter by.

    Returns:
        str: The organism name corresponding to the given organism ID.
    """
    organism = filter_organism_by_id(organism_id)
    if organism is not None:
        return organism.label
    else:
        return ""

def create_orthology_df(orthomap_file: Path):
    """
    Create a DataFrame from the orthologous mapping file.

    Args:
        orthomap_file (Path): The path to the orthologous mapping file in HDF5 format.

    Returns:
        pd.DataFrame: The DataFrame containing the orthologous mapping data.
    """
    try:
        orthology_df = pd.read_hdf(str(orthomap_file))   # do not need to pass in key, as it is the default (only) key
    except Exception as e:
        error_message = f"Error reading orthologous mapping file {orthomap_file}: {e}"
        print(error_message, file=sys.stderr)
        # print stack trace
        import traceback
        traceback.print_exc()
        raise Exception(error_message)
    return orthology_df

def map_dataframe_genes(orig_df: pd.DataFrame, orthomap_file: Path):
    """
    Maps the genes in the original DataFrame to their orthologous genes based on the orthology mapping file.

    Args:
        orig_df (pd.DataFrame): The original DataFrame containing gene data.
        orthomap_file (Path): The path to the orthology mapping file.

    Returns:
        pd.DataFrame: The DataFrame with the genes mapped to their orthologous genes.
    """

    orthology_df = create_orthology_df(orthomap_file)

    # NOTE: Not all genes can be mapped. Unmappable genes do not change in the original DataFrame.

    def get_best_match(gene_id: str) -> str:
        """
        Get the best match for the given gene ID based on the orthology mapping.

        Args:
            gene_id (str): The gene ID to find the best match for.

        Returns:
            str: The best matching gene ID or the original gene ID if no match is found.
        """
        matches = orthology_df[orthology_df["id1"] == gene_id]
        if matches.empty:
            return gene_id

        best_match = matches.sort_values(by="algorithms_match_count", ascending=False).iloc[0]
        return best_match["id2"]

    # Rename orig_df index (id1) using orthologous id (id2).  If id1 maps to multiple id2, then use the id2 with the highest algorithms_match_count
    orig_df.index = orig_df.index.map(get_best_match)

    return orig_df

def map_single_gene(gene_symbol: str, orthomap_file: Path) -> list:
    """
    Maps a single gene symbol to its corresponding orthologous gene symbol(s) using an orthology mapping file.

    Args:
        gene_symbol (str): The gene symbol to be mapped.
        orthomap_file (Path): The path to the orthology mapping file.

    Returns:
        list or None: The list of orthologous gene symbol(s) corresponding to the input gene symbol. Returns an empty list if the gene symbol cannot be mapped.
    """
    orthology_df = create_orthology_df(orthomap_file)

    # Create lowercase gs1 and gs2 columns
    orthology_df["lc_gs1"] = orthology_df["gs1"].str.lower()
    orthology_df["lc_gs2"] = orthology_df["gs2"].str.lower()

    # Check if case-insensitive gene symbol is in dictionary
    lc_gene_symbol = gene_symbol.lower()

    # Return the list of all orthologous gene symbols
    orthologous_genes = orthology_df[orthology_df["lc_gs1"] == lc_gene_symbol]["gs2"].tolist()
    return orthologous_genes

def map_multiple_genes(gene_symbols: list, orthomap_file: Path) -> dict:
    """
    Maps multiple gene symbols to their corresponding orthology gene symbols using an orthology map file.

    Args:
        gene_symbols (list): List of gene symbols to be mapped.
        orthomap_file (Path): Path to the orthology map file.

    Returns:
        dict: A dictionary where the keys are the input gene symbols and the values are lists of orthology gene symbols.
              If a gene symbol does not have any orthology gene symbols, its value is set to an empty list.
    """
    orthology_df = create_orthology_df(orthomap_file)

    # Create lowercase gs1 and gs2 columns
    orthology_df["lc_gs1"] = orthology_df["gs1"].str.lower()
    orthology_df["lc_gs2"] = orthology_df["gs2"].str.lower()

    # Create a dictionary of gene symbol to orthologous gene symbols
    # This checks if case-insensitive gene symbols are in orthology df.
    gene_symbol_dict = {
        gene_symbol: orthology_df[orthology_df["lc_gs1"] == gene_symbol.lower()]["gs2"].tolist()
        for gene_symbol in gene_symbols
    }

    return gene_symbol_dict