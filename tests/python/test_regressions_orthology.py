"""Regression tests for closed GitHub issues in gear.orthology (ortholog lookups and organism helpers)."""

import geardb
import pandas as pd
import pytest
from gear import orthology

ORTHOMAP_COLUMNS = ["id1", "gs1", "id2", "gs2", "algorithms_match_count"]

# Rows shaped like www/feature_mapping/orthomap.*.hdf5: one source gene can map to several
# target genes, and some gs2 values hold several comma-separated symbols.
ORTHOMAP_ROWS = [
    ("ENSDARG01", "calm1a", "ENSG01", "CALM1,CALM3,CALM2", 9),
    ("ENSDARG02", "pou4f3", "ENSG02", "POU4F3", 10),
    ("ENSDARG03", "sox2", "ENSG03", "SOX2", 10),
    ("ENSDARG03", "sox2", "ENSG04", None, 1),      # no target symbol
]


@pytest.fixture
def orthomap_df() -> pd.DataFrame:
    return pd.DataFrame(ORTHOMAP_ROWS, columns=ORTHOMAP_COLUMNS)


@pytest.fixture
def orthomap_file(tmp_path, orthomap_df, monkeypatch):
    """A path whose contents create_orthology_df returns as orthomap_df (no PyTables needed)."""
    path = tmp_path / "orthomap.7.ensembl__1.ensembl.hdf5"
    monkeypatch.setattr(orthology, "create_orthology_df", lambda _path: orthomap_df.copy())
    return path


# ---------------------------------------------------------------- #1150 comma-separated gs2


def test_split_comma_separated_genes_flattens_and_skips_non_strings():
    """Regression for #1150: "CALM1,CALM3,CALM2" splits into three symbols, and None/NaN gs2 values are dropped."""
    genes = orthology.split_comma_separated_genes(["CALM1,CALM3,CALM2", None, float("nan"), "SOX2"])
    assert genes == ["CALM1", "CALM3", "CALM2", "SOX2"]


def test_map_single_gene_splits_comma_separated_orthologs(orthomap_file):
    """Regression for #1150: map_single_gene returns each symbol of a comma-separated gs2 value, never None."""
    assert orthology.map_single_gene("CALM1A", orthomap_file) == ["CALM1", "CALM3", "CALM2"]
    assert orthology.map_single_gene("sox2", orthomap_file) == ["SOX2"]   # the None gs2 row is dropped
    assert orthology.map_single_gene("nope", orthomap_file) == []


def test_map_multiple_genes_splits_comma_separated_orthologs(orthomap_file):
    """Regression for #1150: map_multiple_genes splits comma-separated gs2 values and has no None entries."""
    mapped = orthology.map_multiple_genes(["calm1a", "Pou4f3", "sox2", "nope"], orthomap_file)
    assert mapped == {"calm1a": ["CALM1", "CALM3", "CALM2"], "Pou4f3": ["POU4F3"], "sox2": ["SOX2"], "nope": []}
    assert all(None not in genes for genes in mapped.values())


def test_map_single_gene_reads_real_hdf5_layout(tmp_path, orthomap_df):
    """Regression for #1150: the split also works when reading an HDF5 file laid out like www/feature_mapping (key "orthomap")."""
    pytest.importorskip("tables")
    path = tmp_path / "orthomap.7.ensembl__1.ensembl.hdf5"
    orthomap_df.dropna().to_hdf(path, key="orthomap", mode="w")   # fixed format, as the real files
    assert orthology.map_single_gene("calm1a", path) == ["CALM1", "CALM3", "CALM2"]
    assert orthology.map_multiple_genes(["sox2"], path) == {"sox2": ["SOX2"]}


# ---------------------------------------------------------------- #865 Organism objects


@pytest.fixture
def organisms(monkeypatch):
    orgs = [geardb.Organism(id=1, label="Mouse"), geardb.Organism(id=7, label="Chicken")]
    monkeypatch.setattr(orthology, "organisms", orgs)
    return orgs


def test_filter_organism_by_id_uses_organism_attributes(organisms):
    """Regression for #865: filter_organism_by_id reads Organism.id (not organism["id"]), returning None for an unknown id."""
    assert orthology.filter_organism_by_id(7) is organisms[1]
    assert orthology.filter_organism_by_id(99) is None


def test_get_organism_name_by_id_uses_organism_attributes(organisms):
    """Regression for #865: get_organism_name_by_id reads Organism.label, returning "" for an unknown id."""
    assert orthology.get_organism_name_by_id(1) == "Mouse"
    assert orthology.get_organism_name_by_id(99) == ""
