# orthology.py - functions related to orthology mapping

import sys
import pandas as pd

from pathlib import Path
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
    abs_path_www = Path(__file__).resolve().parents[3] # get www directory
    ORTHOLOG_BASE_DIR = abs_path_www.joinpath("feature_mapping")

    orthomap_file_base = "orthomap.{0}.{2}__{1}.{2}.hdf5".format(gene_organism_id, dataset_organism_id, annotation_source)
    orthomap_file = ORTHOLOG_BASE_DIR.joinpath(orthomap_file_base)
    if not orthomap_file.is_file():
        raise FileNotFoundError("Orthologous mapping file not found: {0}".format(orthomap_file))
    return orthomap_file

def get_ortholog_files_from_dataset(dataset_organism_id: str, annotation_source: str="ensembl"):
    abs_path_www = Path(__file__).resolve().parents[3] # get www directory
    ORTHOLOG_BASE_DIR = abs_path_www.joinpath("feature_mapping")

    orthomap_file_base = "orthomap.*__{0}.{1}.hdf5".format(dataset_organism_id, annotation_source)
    # Find all orthologous mapping files for the given dataset organism ID
    orthomap_files = ORTHOLOG_BASE_DIR.glob(orthomap_file_base)
    if not orthomap_files:
        raise FileNotFoundError("Orthologous mapping files not found for dataset organism {0}".format(get_organism_name_by_id(dataset_organism_id)))

def create_orthology_dict(orthomap_file: Path):
    """
    Create a dictionary of orthologous gene symbols from the orthologous mapping file.

    Args:
        orthomap_file (Path): The path to the orthologous mapping file.

    Returns:
        dict: A dictionary where the keys are gene symbols and the values are the corresponding orthologous gene symbols.
    """
    # Read HDF5 file using Pandas read_hdf
    try:
        orthomap_df = pd.read_hdf(str(orthomap_file))
    except Exception as e:
        raise Exception("Error reading orthologous mapping file: {0}".format(e))
    # Index -> gs1 / id2 / gs2
    orthomap_dict = orthomap_df.to_dict()["id2"]
    # NOTE: Not all genes can be mapped. Unmappable genes do not change in the original dataframe.
    return orthomap_dict

def map_dataframe_genes(orig_df: pd.DataFrame, orthomap_file: Path):
    """
    Remap the passed-in DataFrame to have gene indexes from the orthologous mapping file.

    Parameters:
        orig_df (pd.DataFrame): The original DataFrame to be remapped. The DataFrame's index should be gene symbols.
        orthomap_file (Path): The file containing the orthologous mapping.

    Returns:
        pd.DataFrame: The remapped DataFrame with gene indexes from the orthologous mapping file.
    """
    orthomap_dict = create_orthology_dict(orthomap_file)
    # NOTE: Not all genes can be mapped. Unmappable genes do not change in the original DataFrame.
    return orig_df.rename(index=orthomap_dict)

def map_single_gene(gene_symbol:str, orthomap_file: Path):
    """Map a single gene symbol to the orthologous gene symbol.

    Args:
        gene_symbol (str): The gene symbol to be mapped.
        orthomap_file (Path): The path to the orthology mapping file.

    Returns:
        str: The orthologous gene symbol.

    Raises:
        KeyError: If the gene symbol cannot be mapped.
    """
    # Read HDF5 file using Pandas read_hdf
    orthomap_dict = create_orthology_dict(orthomap_file)
    # NOTE: Not all genes can be mapped. Unmappable genes do not change in the original dataframe.
    return orthomap_dict[gene_symbol]

def map_multiple_genes(gene_symbols:list, orthomap_file: Path):
    """
    Maps multiple gene symbols to their corresponding orthologous gene symbols using an orthology mapping file.

    Args:
        gene_symbols (list): A list of gene symbols to be mapped.
        orthomap_file (Path): The path to the orthology mapping file.

    Returns:
        dict: A dictionary mapping each input gene symbol to its corresponding orthologous gene symbol.
    """
    # Read HDF5 file using Pandas read_hdf
    orthomap_dict = create_orthology_dict(orthomap_file)
    # NOTE: Not all genes can be mapped. Unmappable genes do not change in the original dataframe.
    return { gene_symbol: orthomap_dict[gene_symbol] for gene_symbol in gene_symbols}