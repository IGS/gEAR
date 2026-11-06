"""gear.ontology: ontology term names for sample metadata (#422)."""

import anndata
import pandas as pd
import pytest
import requests
from gear import ontology

LABELS = {
    "UBERON:0002227": "spiral organ of cochlea",
    "UBERON:0000178": "blood",
    "CL:0000601": "cochlear outer hair cell",
    "CL:0011113": "spiral ganglion neuron",
}


@pytest.fixture
def fake_lookup(monkeypatch):
    """No database or OLS: known terms come from LABELS. Records which terms were asked for."""
    asked = []

    def fetch(term_ids):
        asked.extend(term_ids)
        return {t: LABELS[t] for t in term_ids if t in LABELS}

    monkeypatch.setattr(ontology, "_read_cache", lambda term_ids: {})
    monkeypatch.setattr(ontology, "_write_cache", lambda labels: None)
    monkeypatch.setattr(ontology, "_fetch_from_ols", fetch)
    return asked


@pytest.mark.parametrize("value, expected", [
    ("UBERON:0002227", "UBERON:0002227"),
    (" uberon_0002227 ", "UBERON:0002227"),
    ("cl:0000601", "CL:0000601"),
    ("blood", None),
    ("UBERON:", None),
])
def test_normalize_term_id(value, expected):
    assert ontology.normalize_term_id(value) == expected


def test_lookup_terms_uses_cache_first(monkeypatch):
    monkeypatch.setattr(ontology, "_read_cache", lambda term_ids: {"CL:0000601": "cached name"})
    written = []
    monkeypatch.setattr(ontology, "_write_cache", written.append)
    monkeypatch.setattr(ontology, "_fetch_from_ols", lambda term_ids: {t: "from OLS" for t in term_ids if t != "CL:9999999"})

    labels = ontology.lookup_terms(["CL:0000601", "cl_0011113", "CL:9999999", "not a term"])
    assert labels == {"CL:0000601": "cached name", "CL:0011113": "from OLS", "CL:9999999": None}
    assert written == [{"CL:0011113": "from OLS"}]


class FakeResponse:
    def __init__(self, status_code, payload=None):
        self.status_code = status_code
        self._payload = payload or {}

    def json(self):
        return self._payload


def test_fetch_from_ols(monkeypatch):
    requested = []

    def get(self, url, params=None, timeout=None):
        requested.append((url, params["obo_id"]))
        if params["obo_id"] == "CL:0000601":
            return FakeResponse(200, {"_embedded": {"terms": [{"label": "cochlear outer hair cell"}]}})
        return FakeResponse(404)

    monkeypatch.setattr(requests.Session, "get", get)
    assert ontology._fetch_from_ols(["CL:0000601", "CL:9999999"]) == {"CL:0000601": "cochlear outer hair cell"}
    assert requested[0] == (f"{ontology.OLS_API_URL}/ontologies/cl/terms", "CL:0000601")


def test_fetch_from_ols_stops_when_unreachable(monkeypatch):
    calls = []

    def get(self, url, params=None, timeout=None):
        calls.append(params["obo_id"])
        raise requests.ConnectionError("no network")

    monkeypatch.setattr(requests.Session, "get", get)
    assert ontology._fetch_from_ols(["CL:0000601", "CL:0011113"]) == {}
    assert calls == ["CL:0000601"]


def test_names_are_added_for_term_ids(fake_lookup):
    obs = pd.DataFrame({
        "tissue_ontology": ["UBERON:0002227", "UBERON:0002227", "UBERON:0000178"],
        "cell_ontology": ["CL:0000601", "CL:0000601, CL:0011113", None],
    }, index=["c1", "c2", "c3"])

    obs, summary = ontology.annotate_sample_metadata(obs, recommended_columns=[])

    assert list(obs["tissue_ontology_name"]) == ["spiral organ of cochlea", "spiral organ of cochlea", "blood"]
    assert list(obs["cell_ontology_name"]) == ["cochlear outer hair cell",
                                               "cochlear outer hair cell, spiral ganglion neuron", ""]
    assert isinstance(obs["cell_ontology_name"].dtype, pd.CategoricalDtype)
    assert summary == {
        "ontology_terms": {"tissue_type": ["UBERON:0000178", "UBERON:0002227"],
                           "cell_type": ["CL:0000601", "CL:0011113"]},
        "warnings": [],
    }
    # Each distinct term is looked up once
    assert sorted(fake_lookup) == sorted(set(fake_lookup))


def test_column_aliases_are_recognized(fake_lookup):
    obs = pd.DataFrame({"Cell_Type_Ontology_Term_ID": ["CL:0000601"], "UBERON_ID": ["UBERON:0000178"]})
    obs, summary = ontology.annotate_sample_metadata(obs, recommended_columns=[])
    assert summary["ontology_terms"] == {"tissue_type": ["UBERON:0000178"], "cell_type": ["CL:0000601"]}
    assert obs["cell_ontology_name"].tolist() == ["cochlear outer hair cell"]


def test_bad_and_unknown_terms_are_warnings(fake_lookup):
    obs = pd.DataFrame({"tissue_ontology": ["blood", "CL:0000601", "UBERON:9999999", "UBERON:0000178"]})
    obs, summary = ontology.annotate_sample_metadata(obs, recommended_columns=[])

    assert summary["ontology_terms"] == {"tissue_type": ["UBERON:0000178", "UBERON:9999999"]}
    assert len(summary["warnings"]) == 3
    assert "blood" in summary["warnings"][0]
    assert "CL:0000601" in summary["warnings"][1]
    assert "UBERON:9999999" in summary["warnings"][2]
    # Unresolved values keep their ID as the name
    assert obs["tissue_ontology_name"].tolist() == ["blood", "CL:0000601", "UBERON:9999999", "blood"]


def test_missing_recommended_columns(fake_lookup):
    obs = pd.DataFrame({"Treatment": ["ctl"], "tissue_ontology_term_id": ["UBERON:0000178"]})
    _, summary = ontology.annotate_sample_metadata(
        obs, recommended_columns=["tissue_ontology", "treatment", "age_value", "sex_assigned_at_birth"])
    assert summary["warnings"] == [
        "These recommended sample metadata columns are missing: age_value, sex_assigned_at_birth. "
        "See the sample metadata template for details."
    ]


@pytest.mark.parametrize("columns, warned", [
    (["UMAP_1", "UMAP_2"], False),
    (["PC_1", "PC_2"], False),
    (["cell_clusters"], False),
    (["seurat_clusters"], False),
    (["RNA_snn_res.0.8"], False),
    # One coordinate column alone isn't an embedding gEAR can plot
    (["UMAP_1"], True),
    (["condition"], True),
])
def test_primary_analysis_columns_checked(fake_lookup, monkeypatch, columns, warned):
    monkeypatch.setattr(ontology, "get_recommended_obs_columns", lambda: [])
    metadata = {"primary_analysis_completed": "Yes"}
    ontology.annotate_upload_obs(anndata.AnnData(obs=pd.DataFrame({c: [1] for c in columns})), metadata)
    assert bool(metadata["sample_metadata_warnings"]) is warned


def test_primary_analysis_not_checked_when_not_claimed(fake_lookup, monkeypatch):
    monkeypatch.setattr(ontology, "get_recommended_obs_columns", lambda: [])
    metadata = {"primary_analysis_completed": "No"}
    ontology.annotate_upload_obs(anndata.AnnData(obs=pd.DataFrame({"condition": [1]})), metadata)
    assert metadata["sample_metadata_warnings"] == []


def test_recommended_columns_from_config(monkeypatch):
    class FakeConfig:
        def __init__(self, value):
            self.value = value

        def parse(self):
            return self

        def get(self, section, option, fallback=None):
            return self.value if self.value is not None else fallback

    monkeypatch.setattr(ontology, "ServerConfig", lambda: FakeConfig("treatment, age_value"))
    assert ontology.get_recommended_obs_columns() == ["treatment", "age_value"]
    monkeypatch.setattr(ontology, "ServerConfig", lambda: FakeConfig(""))
    assert ontology.get_recommended_obs_columns() == []
    monkeypatch.setattr(ontology, "ServerConfig", lambda: FakeConfig(None))
    assert ontology.get_recommended_obs_columns() == list(ontology.DEFAULT_RECOMMENDED_OBS_COLUMNS)


def test_annotate_upload_obs_records_results(fake_lookup, monkeypatch):
    monkeypatch.setattr(ontology, "get_recommended_obs_columns", lambda: ["age_value"])
    metadata = {"primary_analysis_completed": "No"}
    adata = anndata.AnnData(obs=pd.DataFrame({"cell_ontology": ["CL:0000601"]}))
    ontology.annotate_upload_obs(adata, metadata)
    assert "cell_ontology_name" in adata.obs
    assert metadata["ontology_terms"] == {"cell_type": ["CL:0000601"]}
    assert len(metadata["sample_metadata_warnings"]) == 1


def test_annotate_upload_obs_never_raises(monkeypatch):
    def broken(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(ontology, "annotate_sample_metadata", broken)
    monkeypatch.setattr(ontology, "detect_primary_analysis", broken)
    metadata = {"primary_analysis_completed": "Yes"}
    obs = pd.DataFrame({"cell_ontology": ["CL:0000601"]})
    adata = anndata.AnnData(obs=obs)
    ontology.annotate_upload_obs(adata, metadata)
    assert list(adata.obs.columns) == ["cell_ontology"]
    assert metadata["ontology_terms"] == {} and len(metadata["sample_metadata_warnings"]) == 1


def test_find_standard_obs_columns():
    obs = pd.DataFrame({
        "Tissue_Ontology_ID": ["UBERON:0002227"], "SEX_ASSIGNED_AT_BIRTH": ["female"], "treatment": ["ctr"],
        "seurat_clusters": ["0"], "UMAP_1": [1.0], "UMAP_2": [2.0], "pca_1": [0.1], "pca_2": [0.2],
        "orig.ident": ["sample1"], "batch": ["a"],
    })
    assert ontology.find_standard_obs_columns(obs) == {
        "tissue_ontology": "Tissue_Ontology_ID",
        "treatment": "treatment",
        "sex_assigned_at_birth": "SEX_ASSIGNED_AT_BIRTH",
        "cell_clusters": "seurat_clusters",
        "umap": ["UMAP_1", "UMAP_2"],
        "pca": ["pca_1", "pca_2"],
    }
    assert ontology.find_standard_obs_columns(pd.DataFrame({"batch": ["a"]})) == {}


def test_annotate_upload_obs_saves_column_map(fake_lookup, monkeypatch, tmp_path):
    """The column map lands in adata.uns, names columns gEAR added, and survives writing the H5AD."""
    monkeypatch.setattr(ontology, "get_recommended_obs_columns", lambda: [])
    obs = pd.DataFrame({"cell_type_ontology_term_id": ["CL:0000601"], "leiden": ["1"],
                        "tSNE_1": [1.0], "tSNE_2": [2.0]}, index=["c1"])
    adata = anndata.AnnData(obs=obs)
    ontology.annotate_upload_obs(adata, {"primary_analysis_completed": "Yes"})

    expected = {"cell_ontology": "cell_type_ontology_term_id", "cell_ontology_name": "cell_ontology_name",
                "cell_clusters": "leiden", "tsne": ["tSNE_1", "tSNE_2"]}
    assert adata.uns[ontology.OBS_COLUMN_MAP_KEY] == expected

    path = tmp_path / "test.h5ad"
    adata.write_h5ad(path)
    saved = anndata.read_h5ad(path).uns[ontology.OBS_COLUMN_MAP_KEY]
    assert {key: value if isinstance(value, str) else list(value) for key, value in saved.items()} == expected
