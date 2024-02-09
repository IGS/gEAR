# orthology.py - functions related to orthology mapping

import sys
import pandas as pd

from pathlib import Path
from collections import defaultdict

# append parent directory to path
sys.path.append(str(Path(__file__).resolve().parents[1]))

import geardb

organisms = geardb.OrganismCollection().get_all()

def filter_organism_by_id(organism_id: str):
    """Filter the organisms list by the given organism ID.

    Args:
        organism_id (str): The organism ID to filter by.

    Returns:
        dict: The organism dictionary corresponding to the given organism ID.
    """
    return next((item for item in organisms if item["id"] == organism_id), None)

def get_organism_name_by_id(organism_id: str):
    """Get the organism name corresponding to the given organism ID.

    Args:
        organism_id (str): The organism ID to filter by.

    Returns:
        str: The organism name corresponding to the given organism ID.
    """
    return filter_organism_by_id(organism_id)["name"]

def get_ortholog_file(gene_organism_id: str, dataset_organism_id: str, annotation_source: str="ensembl"):
    """
    Get the path to the orthologous mapping file based on the given gene organism ID, dataset organism ID, and annotation source.

    Args:
        gene_organism_id (str): The ID of the gene organism.
        dataset_organism_id (str): The ID of the dataset organism.
        annotation_source (str, optional): The source of the annotation. Defaults to "ensembl".

    Returns:
        Path: The path to the orthologous mapping file.

    Raises:
        FileNotFoundError: If the orthologous mapping file is not found.
    """
    abs_path_www = Path(__file__).resolve().parents[2].joinpath("www") # get www directory
    ORTHOLOG_BASE_DIR = abs_path_www.joinpath("feature_mapping")

    orthomap_file_base = "orthomap.{0}.{2}__{1}.{2}.hdf5".format(dataset_organism_id, gene_organism_id, annotation_source)
    orthomap_file = ORTHOLOG_BASE_DIR.joinpath(orthomap_file_base)
    if not orthomap_file.is_file():
        raise FileNotFoundError("Orthologous mapping file not found: {0}".format(orthomap_file))
    return orthomap_file

def get_ortholog_files_from_dataset(dataset_organism_id: str, annotation_source: str="ensembl"):
    """
    Retrieves the orthologous mapping files for a given dataset organism ID.

    Args:
        dataset_organism_id (str): The ID of the dataset organism.
        annotation_source (str, optional): The annotation source. Defaults to "ensembl".

    Returns:
        List[Path]: A list of Path objects representing the orthologous mapping files.

    Raises:
        FileNotFoundError: If no orthologous mapping files are found for the given dataset organism ID.
    """
    abs_path_www = Path(__file__).resolve().parents[2].joinpath("www") # get www directory
    ORTHOLOG_BASE_DIR = abs_path_www.joinpath("feature_mapping")

    # Find all orthologous mapping files for the given dataset organism ID
    orthomap_file_base = "orthomap.{0}.{1}__[0-9].{1}.hdf5".format(dataset_organism_id, annotation_source)
    orthomap_files = list(ORTHOLOG_BASE_DIR.glob(orthomap_file_base))

    # sort by gene organism id
    orthomap_files.sort(key=lambda x: int(x.name.split("__")[1].split(".")[0]))
    if not orthomap_files:
        raise FileNotFoundError("Orthologous mapping files not found for dataset organism {0}".format(get_organism_name_by_id(dataset_organism_id)))
    return orthomap_files


def create_orthology_dict(orthomap_file: Path):
    """
    Create a dictionary mapping gene IDs to a list of corresponding orthologous gene symbols.

    Args:
        orthomap_file (Path): The path to the HDF5 file containing the orthologous mapping data.

    Returns:
        dict: A dictionary where the keys are Ensembl IDs (id2) from the searched organism and the values are lists of corresponding orthologous Ensembl IDs (id1) from the dataset organism.
    """

    # Read HDF5 file using Pandas read_hdf
    try:
        orthomap_df = pd.read_hdf(str(orthomap_file))
    except Exception as e:
        raise Exception("Error reading orthologous mapping file: {0}".format(e))
    # Index (id1) -> gs1 / id2 / gs2
    # NOTE: Not all genes can be mapped. Unmappable genes have an empty list

    # for each index, get list of id1 values
    orthomap_dict = defaultdict(list)
    for index, row in orthomap_df.iterrows():
        orthomap_dict[row["id2"]].append(row["id1"])
    return orthomap_dict

def create_orthology_gene_symbol_dict(orthomap_file: Path):
    """
    Create a dictionary mapping orthology gene symbols.

    Args:
        orthomap_file (Path): The path to the orthologous mapping file.

    Returns:
        dict: A dictionary where the keys are gene symbols (gs2) from the searched organism and the values are lists of corresponding orthologous gene symbols (gs1) from the dataset organism.
    """
    # Read HDF5 file using Pandas read_hdf
    try:
        orthomap_df = pd.read_hdf(str(orthomap_file))
    except Exception as e:
        raise Exception("Error reading orthologous mapping file: {0}".format(e))
    # Index (id1) -> gs1 / id2 / gs2
    # NOTE: Not all genes can be mapped. Unmappable genes have an empty list

    # for each index, get list of gs2 values
    # gs2 is the gene we searched for so need to treat gs1 as isoforms
    orthomap_dict = defaultdict(list)
    for index, row in orthomap_df.iterrows():
        orthomap_dict[row["gs2"]].append(row["gs1"])
    return orthomap_dict

def map_dataframe_genes(orig_df: pd.DataFrame, orthomap_file: Path):
    """
    Maps the gene symbols in the original DataFrame to their orthologous gene symbols using the provided orthology mapping file.

    Args:
        orig_df (pd.DataFrame): The original DataFrame containing gene symbols.
        orthomap_file (Path): The path to the orthology mapping file.

    Returns:
        pd.DataFrame: The DataFrame with gene symbols mapped to their orthologous gene symbols.
    """
    orthomap_dict = create_orthology_dict(orthomap_file)
    # NOTE: Not all genes can be mapped. Unmappable genes do not change in the original DataFrame.

    # Rename index using orthologous gene symbols.  If more than one orthologous gene symbol is found, the first one is used.
    # This is primarily used for projections, but we could take all orthologs as an option if a scenario arises
    for index in orig_df.index:
        mapped_gene_symbols = orthomap_dict.get(index, None)
        orig_df.rename(index={index: mapped_gene_symbols}, inplace=True)
    return orig_df

def map_single_gene(gene_symbol:str, orthomap_file: Path):
    """
    Maps a single gene symbol to its corresponding orthologous gene symbol(s) using an orthology mapping file.

    Args:
        gene_symbol (str): The gene symbol to be mapped.
        orthomap_file (Path): The path to the orthology mapping file.

    Returns:
        list or None: The list of orthologous gene symbol(s) corresponding to the input gene symbol. Returns an empty list if the gene symbol cannot be mapped.
    """
    # Read HDF5 file using Pandas read_hdf
    gene_symbol_dict = create_orthology_gene_symbol_dict(orthomap_file)

    # Convert dictionary keys to lowercase
    gene_symbol_dict = {key.lower(): value for key, value in gene_symbol_dict.items()}

    # Check if case-insensitive gene symbol is in dictionary
    lc_gene_symbol = gene_symbol.lower()

    # returns the list of all isoforms
    return gene_symbol_dict.get(lc_gene_symbol, [])

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
    # Read HDF5 file using Pandas read_hdf
    gene_symbol_dict = create_orthology_gene_symbol_dict(orthomap_file)

    # Convert dictionary keys to lowercase
    gene_symbol_dict = {key.lower(): value for key, value in gene_symbol_dict.items()}

    # Check if case-insensitive gene symbols are in dictionary.
    case_insensitive = {gene_symbol: gene_symbol_dict.get(gene_symbol.lower(), []) for gene_symbol in gene_symbols}

    # returns the list of all isoforms for each original gene symbol
    return case_insensitive