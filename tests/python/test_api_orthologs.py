"""POST /h5ad/<dataset_id>/orthologs and its mapping helpers in resources.orthologs."""

import types

import anndata
import numpy as np
import pandas as pd
import pytest
from flask import Flask
from flask_restful import Api
from resources import orthologs

DATASET_GENES = {"Sox2", "Pou4f3", "Calm1"}


@pytest.fixture
def gene_map():
    data = orthologs.build_gene_map(DATASET_GENES)
    return data["gene_map_set"], data["normalize_gene"]


# ---------------------------------------------------------------- pure helpers


def test_build_gene_map_is_case_insensitive(gene_map):
    gene_set, normalize = gene_map
    assert gene_set == {"sox2", "pou4f3", "calm1"}
    assert normalize("SOX2") == "Sox2"
    assert normalize("nope") is None


def test_check_gene_in_dataset_lowercases_input():
    assert orthologs.check_gene_in_dataset({"tp53", "brca1"}, "TP53")
    assert orthologs.check_gene_in_dataset({"123"}, 123)
    assert not orthologs.check_gene_in_dataset({"tp53"}, "brca1")


def test_normalize_mapped_genes_keeps_dataset_genes_in_dataset_case(gene_map):
    gene_set, normalize = gene_map
    assert orthologs.normalize_mapped_genes(["SOX2", "NOPE", "calm1"], gene_set, normalize) == ["Sox2", "Calm1"]
    assert orthologs.normalize_mapped_genes([], gene_set, normalize) == []


def test_last_chance_map(gene_map):
    gene_set, normalize = gene_map
    assert orthologs.last_chance_map("POU4F3", gene_set, normalize) == ["Pou4f3"]
    assert orthologs.last_chance_map("nope", gene_set, normalize) == []


# ---------------------------------------------------------------- ortholog files


@pytest.fixture
def ortholog_files(monkeypatch):
    """Ortholog file lookups keyed by organism; records each lookup."""
    calls = []
    gene_files = {(7, 1): "orthomap.7.ensembl__1.ensembl.hdf5"}
    dataset_files = {1: ["orthomap.2.ensembl__1.ensembl.hdf5", "orthomap.3.ensembl__1.ensembl.hdf5"]}

    def get_ortholog_file(gene_org, dataset_org, annotation):
        calls.append(("gene", gene_org, dataset_org))
        if (gene_org, dataset_org) not in gene_files:
            raise FileNotFoundError(f"no file for {gene_org} -> {dataset_org}")
        return gene_files[(gene_org, dataset_org)]

    def get_ortholog_files_from_dataset(dataset_org, annotation):
        calls.append(("dataset", dataset_org))
        if dataset_org not in dataset_files:
            raise FileNotFoundError(f"no files for {dataset_org}")
        return list(dataset_files[dataset_org])

    monkeypatch.setattr(orthologs, "get_ortholog_file", get_ortholog_file)
    monkeypatch.setattr(orthologs, "get_ortholog_files_from_dataset", get_ortholog_files_from_dataset)
    return calls


def test_fetch_ortholog_files_same_organism_uses_dataset_files(ortholog_files):
    files = orthologs.fetch_ortholog_files(1, 1)
    assert files == ["orthomap.2.ensembl__1.ensembl.hdf5", "orthomap.3.ensembl__1.ensembl.hdf5"]
    assert ortholog_files == [("dataset", 1)]
    assert orthologs.fetch_ortholog_files(None, 1) == files


def test_fetch_ortholog_files_cross_organism_puts_gene_file_first(ortholog_files):
    files = orthologs.fetch_ortholog_files(7, 1)
    assert files == ["orthomap.7.ensembl__1.ensembl.hdf5",
                     "orthomap.2.ensembl__1.ensembl.hdf5", "orthomap.3.ensembl__1.ensembl.hdf5"]


def test_fetch_ortholog_files_exclusive_only_uses_gene_file(ortholog_files):
    assert orthologs.fetch_ortholog_files(7, 1, exclusive_org=True) == ["orthomap.7.ensembl__1.ensembl.hdf5"]
    assert ("dataset", 1) not in ortholog_files


def test_fetch_ortholog_files_reraises_missing_file(ortholog_files):
    with pytest.raises(FileNotFoundError):
        orthologs.fetch_ortholog_files(9, 1)
    with pytest.raises(FileNotFoundError):
        orthologs.fetch_ortholog_files(None, 99)


def test_get_mapped_gene_symbols_uses_first_file_with_a_mapping(ortholog_files, monkeypatch):
    mappings = {
        "orthomap.2.ensembl__1.ensembl.hdf5": {"sox2": []},
        "orthomap.3.ensembl__1.ensembl.hdf5": {"sox2": ["Sox2"]},
    }
    monkeypatch.setattr(orthologs, "map_multiple_genes", lambda genes, path: mappings[path])
    assert orthologs.get_mapped_gene_symbols(["sox2"], None, 1) == {"sox2": ["Sox2"]}


def test_get_mapped_gene_symbols_without_files_or_mapping(ortholog_files, monkeypatch):
    monkeypatch.setattr(orthologs, "map_multiple_genes", lambda genes, path: {g: [] for g in genes})
    assert orthologs.get_mapped_gene_symbols(["sox2"], None, 99) == {"sox2": []}   # no files
    assert orthologs.get_mapped_gene_symbols(["sox2"], None, 1) == {"sox2": []}
    with pytest.raises(Exception, match="No orthologous mapping"):
        orthologs.get_mapped_gene_symbols(["sox2"], 7, 1, exclusive_org=True)


# ---------------------------------------------------------------- map_all_genes


@pytest.fixture
def ortho_lookups(monkeypatch):
    """get_mapped_gene_symbols returns ortho_lookups.result (or raises it); records the genes it was asked for."""
    lookups = types.SimpleNamespace(calls=[], result={})

    def fake(gene_symbols, gene_organism_id, dataset_organism_id, exclusive_org=False):
        lookups.calls.append(list(gene_symbols))
        if isinstance(lookups.result, Exception):
            raise lookups.result
        return lookups.result

    monkeypatch.setattr(orthologs, "get_mapped_gene_symbols", fake)
    return lookups


def map_all(genes, gene_map, gene_org=1, dataset_org=1, exclusive_org=False):
    gene_set, normalize = gene_map
    return orthologs.map_all_genes(genes, gene_org, dataset_org, gene_set, exclusive_org, normalize)


def test_map_all_genes_same_organism_maps_directly(gene_map, ortho_lookups):
    assert map_all(["SOX2", "calm1"], gene_map) == {"SOX2": ["Sox2"], "calm1": ["Calm1"]}
    assert ortho_lookups.calls == []


def test_map_all_genes_orthology_only_for_unmapped(gene_map, ortho_lookups):
    ortho_lookups.result = {"pou4f3b": ["POU4F3"]}
    assert map_all(["sox2", "pou4f3b"], gene_map) == {"sox2": ["Sox2"], "pou4f3b": ["Pou4f3"]}
    assert ortho_lookups.calls == [["pou4f3b"]]


def test_map_all_genes_cross_organism_falls_back_to_symbol_in_dataset(gene_map, ortho_lookups):
    ortho_lookups.result = {"sox2": ["NOT_IN_DATASET"], "calm1a": ["CALM1"], "zzz": []}
    result = map_all(["sox2", "calm1a", "zzz"], gene_map, gene_org=7)
    assert result == {"sox2": ["Sox2"], "calm1a": ["Calm1"], "zzz": []}
    assert ortho_lookups.calls == [["sox2", "calm1a", "zzz"]]


def test_map_all_genes_orthology_error(gene_map, ortho_lookups):
    ortho_lookups.result = RuntimeError("No orthologous mapping found")
    assert map_all(["zzz"], gene_map, gene_org=7) == {"zzz": []}
    with pytest.raises(RuntimeError):
        map_all(["zzz"], gene_map, gene_org=7, exclusive_org=True)


# ---------------------------------------------------------------- POST endpoint


def genes_adata():
    var = pd.DataFrame({"gene_symbol": sorted(DATASET_GENES)}, index=["E1", "E2", "E3"])
    return anndata.AnnData(X=np.zeros((2, 3), dtype=np.float32), var=var)


@pytest.fixture
def client(monkeypatch, ortho_lookups):
    datasets = {
        "DS1": types.SimpleNamespace(id="DS1", dtype="single-cell-rnaseq", organism_id=1),
        "SP1": types.SimpleNamespace(id="SP1", dtype="spatial", organism_id=1),
        "NOORG": types.SimpleNamespace(id="NOORG", dtype="single-cell-rnaseq", organism_id=None),
    }
    loads = []

    def loader(name):
        def _load(analysis_id, dataset_id, session_id):
            loads.append((name, analysis_id, dataset_id, session_id))
            if isinstance(test_client.adata, Exception):
                raise test_client.adata
            return test_client.adata
        return _load

    monkeypatch.setattr(orthologs.geardb, "get_dataset_by_id", lambda dataset_id: datasets.get(dataset_id))
    monkeypatch.setattr(orthologs, "get_adata_shadow", loader("get_adata_shadow"))
    monkeypatch.setattr(orthologs, "get_spatial_adata", loader("get_spatial_adata"))
    exclusive = []
    real_map_all = orthologs.map_all_genes

    def spy_map_all(*args):
        exclusive.append(args[4])
        return real_map_all(*args)

    monkeypatch.setattr(orthologs, "map_all_genes", spy_map_all)

    app = Flask(__name__)
    Api(app).add_resource(orthologs.Orthologs, "/h5ad/<dataset_id>/orthologs")
    test_client = app.test_client()
    test_client.adata = genes_adata()
    test_client.loads = loads
    test_client.exclusive = exclusive
    return test_client


def post(client, body, dataset_id="DS1"):
    return client.post(f"/h5ad/{dataset_id}/orthologs", json=body)


def test_post_maps_genes(client, ortho_lookups):
    ortho_lookups.result = {"calm1a": ["CALM1"]}
    client.set_cookie("gear_session_id", "sess")
    resp = post(client, {"gene_symbols": ["SOX2", "calm1a"], "gene_organism_id": "1", "analysis": {"id": "A1"}})
    assert resp.status_code == 200
    assert resp.get_json() == {"success": 1, "mapping": {"SOX2": ["Sox2"], "calm1a": ["Calm1"]}}
    assert client.loads == [("get_adata_shadow", "A1", "DS1", "sess")]
    assert client.exclusive == [False]


def test_post_spatial_uses_spatial_loader(client):
    assert post(client, {"gene_symbols": ["sox2"], "gene_organism_id": 1}, dataset_id="SP1").status_code == 200
    assert [load[0] for load in client.loads] == ["get_spatial_adata"]


@pytest.mark.parametrize("value, expected", [
    ("true", True), ("TRUE", True), ("false", False), (True, True), (False, False), (None, False)])
def test_post_exclusive_org_accepts_strings_and_booleans(client, value, expected):
    resp = post(client, {"gene_symbols": ["sox2"], "gene_organism_id": 1, "exclusive_org": value})
    assert resp.status_code == 200, resp.get_data(as_text=True)
    assert client.exclusive == [expected]


@pytest.mark.parametrize("body, dataset_id, error", [
    ({"gene_symbols": []}, "DS1", "No gene symbols provided."),
    ({"gene_symbols": ["sox2"]}, "NOPE", "The dataset was not found."),
    ({"gene_symbols": ["sox2"]}, "NOORG", "The dataset does not have an associated organism."),
])
def test_post_request_errors(client, body, dataset_id, error):
    resp = post(client, body, dataset_id=dataset_id)
    assert resp.status_code == 400
    assert resp.get_json() == {"error": error}
    assert client.loads == []


def test_post_empty_json_body(client):
    resp = client.post("/h5ad/DS1/orthologs", data="{}", content_type="application/json")
    assert resp.status_code == 400
    assert resp.get_json() == {"error": "Invalid JSON body."}


@pytest.mark.parametrize("adata, message", [
    (FileNotFoundError("gone"), "No dataset file found."),
    (anndata.AnnData(X=np.zeros((1, 1), dtype=np.float32)), "Dataset does not contain 'gene_symbol' field."),
])
def test_post_dataset_load_errors(client, adata, message):
    client.adata = adata
    resp = post(client, {"gene_symbols": ["sox2"]})
    assert resp.status_code == 400
    assert resp.get_json() == {"success": -1, "message": message}
