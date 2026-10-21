"""
Small Flask API resources: analyses, availableAnalysisTools, genes, h5ad obs info, HiGlass gene lookup,
SVG expression data and dataset displays.

Resources are called through Flask's test client. Database lookups and AnnData loaders are monkeypatched,
and requests to HiGlass are mocked.
"""

import types
from urllib.parse import parse_qs, urlsplit

import anndata
import numpy as np
import pandas as pd
import pytest
import requests
from flask import Flask
from flask_restful import Api
from gear.plotting import PlotError
from resources import (
    analyses,
    available_analysis_tools,
    dataset_display,
    gene_symbols,
    h5ad,
    higlass,
    svg_data,
)

SINGLE_CELL = types.SimpleNamespace(id="DS1", dtype="single-cell-rnaseq")
SPATIAL = types.SimpleNamespace(id="DS1", dtype="spatial")


def make_client(resource, route):
    app = Flask(__name__)
    Api(app).add_resource(resource, route)
    return app.test_client()


def ok_json(resp):
    assert resp.status_code == 200, resp.get_data(as_text=True)
    return resp.get_json()


def record_loaders(monkeypatch, mod, adata=None, error=None):
    """Patch get_adata_shadow/get_spatial_adata on mod to return adata (or raise error); returns the call log."""
    calls = []

    def loader(name):
        def _load(analysis_id, dataset_id, session_id, *args, **kwargs):
            calls.append((name, analysis_id, dataset_id, session_id))
            if error is not None:
                raise error
            return adata
        return _load

    for name in ("get_adata_shadow", "get_spatial_adata"):
        monkeypatch.setattr(mod, name, loader(name))
    return calls


# ---------------------------------------------------------------- analyses


class FakeAnalysis:
    def __init__(self, id, tsne=0, umap=0):
        self.id = id
        self.tsne = {"tsne_calculated": tsne, "umap_calculated": umap}

    def _serialize_json(self):
        return {"id": self.id}


def test_tsne_or_umap_present():
    assert analyses.tsne_or_umap_present(FakeAnalysis("a", tsne=1))
    assert analyses.tsne_or_umap_present(FakeAnalysis("a", umap=1))
    assert not analyses.tsne_or_umap_present(FakeAnalysis("a"))


@pytest.fixture
def analyses_client(monkeypatch):
    lookups = []

    class FakeCollection:
        def get_all_by_dataset_id(self, **kwargs):
            lookups.append(kwargs)
            self.public = [FakeAnalysis("pub-tsne", tsne=1), FakeAnalysis("pub-none")]
            self.user_saved = [FakeAnalysis("mine-umap", umap=1)] if kwargs["user_id"] else []

    users = {"owner": types.SimpleNamespace(id=5)}
    datasets = {"DS1": SINGLE_CELL, "SP1": SPATIAL}
    monkeypatch.setattr(analyses.geardb, "get_user_from_session_id", lambda sid: users.get(sid))
    monkeypatch.setattr(analyses.geardb, "get_dataset_by_id", lambda dataset_id: datasets.get(dataset_id))
    monkeypatch.setattr(analyses, "AnalysisCollection", FakeCollection)
    client = make_client(analyses.Analyses, "/h5ad/<dataset_id>/analyses")
    client.lookups = lookups
    return client


def test_analyses_lists_only_analyses_with_tsne_or_umap(analyses_client):
    analyses_client.set_cookie("gear_session_id", "owner")
    result = ok_json(analyses_client.get("/h5ad/DS1/analyses"))
    assert result == {"success": 1, "public": [{"id": "pub-tsne"}], "private": [{"id": "mine-umap"}]}
    assert analyses_client.lookups == [
        {"user_id": 5, "session_id": "owner", "dataset_id": "DS1", "is_spatial": False}]


def test_analyses_logged_out_has_no_private(analyses_client):
    result = ok_json(analyses_client.get("/h5ad/SP1/analyses"))
    assert result["private"] == []
    assert analyses_client.lookups[0]["user_id"] is None
    assert analyses_client.lookups[0]["is_spatial"] is True


def test_analyses_unknown_dataset_is_error_with_200(analyses_client):
    result = ok_json(analyses_client.get("/h5ad/NOPE/analyses"))
    assert result == {"success": -1, "message": "No dataset found with that ID"}
    assert analyses_client.lookups == []


# ---------------------------------------------------------------- availableAnalysisTools


@pytest.fixture
def tools_client(monkeypatch):
    dtypes = {"S-SC": "single-cell-rnaseq", "S-BULK": "bulk-rnaseq", "S-SVG": "svg-expression", "S-ODD": "epiviz"}
    monkeypatch.setattr(available_analysis_tools.geardb, "get_dtype_by_share_id", lambda share_id: dtypes.get(share_id))
    return make_client(available_analysis_tools.AvailableAnalysisTools, "/h5ad/<share_uid>/availableAnalysisTools")


@pytest.mark.parametrize("share_id, expected", [
    ("S-SC", {"dataset-curator": True, "multigene-viewer": True, "compare-tool": True, "sc-workbench": True}),
    ("S-BULK", {"dataset-curator": True, "multigene-viewer": True, "compare-tool": True, "sc-workbench": False}),
    ("S-SVG", {"dataset-curator": True, "multigene-viewer": False, "compare-tool": False, "sc-workbench": False}),
    ("S-ODD", {"dataset-curator": False, "multigene-viewer": False, "compare-tool": False, "sc-workbench": False}),
])
def test_available_analysis_tools_by_dtype(tools_client, share_id, expected):
    result = ok_json(tools_client.get(f"/h5ad/{share_id}/availableAnalysisTools"))
    assert result == {"success": 1, "share_id": share_id, "available_analysis_tools": expected}


def test_available_analysis_tools_unknown_share_id(tools_client):
    result = ok_json(tools_client.get("/h5ad/NOPE/availableAnalysisTools"))
    assert result == {"success": -1, "message": "No dataset found with that share ID"}


# ---------------------------------------------------------------- genes


def genes_adata(symbols):
    var = pd.DataFrame({"gene_symbol": symbols}, index=[f"E{i}" for i in range(len(symbols))])
    return anndata.AnnData(X=np.zeros((2, len(symbols)), dtype=np.float32), var=var)


@pytest.fixture
def genes_client(monkeypatch):
    datasets = {"DS1": SINGLE_CELL, "SP1": SPATIAL}
    monkeypatch.setattr(gene_symbols.geardb, "get_dataset_by_id", lambda dataset_id: datasets.get(dataset_id))
    return make_client(gene_symbols.GeneSymbols, "/h5ad/<dataset_id>/genes")


def test_gene_symbols_lists_var_gene_symbols(genes_client, monkeypatch):
    calls = record_loaders(monkeypatch, gene_symbols, genes_adata(["Sox2", "Pou4f3"]))
    genes_client.set_cookie("gear_session_id", "sess")
    result = ok_json(genes_client.get("/h5ad/DS1/genes?analysis=A1"))
    assert result == {"success": 1, "gene_symbols": ["Sox2", "Pou4f3"]}
    assert calls == [("get_adata_shadow", "A1", "DS1", "sess")]


def test_gene_symbols_spatial_uses_spatial_loader(genes_client, monkeypatch):
    calls = record_loaders(monkeypatch, gene_symbols, genes_adata(["Sox2"]))
    ok_json(genes_client.get("/h5ad/SP1/genes"))
    assert calls == [("get_spatial_adata", None, "SP1", None)]


@pytest.mark.parametrize("symbols", [
    pd.Categorical(["Sox2", np.nan, "Pou4f3"]),
    pd.array(["Sox2", None, "Pou4f3"], dtype=object),
])
def test_gene_symbols_missing_symbols_become_none(genes_client, monkeypatch, symbols):
    record_loaders(monkeypatch, gene_symbols, genes_adata(symbols))
    result = ok_json(genes_client.get("/h5ad/DS1/genes"))
    assert result["gene_symbols"] == ["Sox2", "None", "Pou4f3"]


@pytest.mark.parametrize("error, message", [
    (FileNotFoundError("gone"), "No file found for this dataset"),
    (ValueError("bad analysis"), "bad analysis"),
])
def test_gene_symbols_loader_errors(genes_client, monkeypatch, error, message):
    record_loaders(monkeypatch, gene_symbols, error=error)
    assert ok_json(genes_client.get("/h5ad/DS1/genes")) == {"success": -1, "message": message}


def test_gene_symbols_unknown_dataset(genes_client):
    assert ok_json(genes_client.get("/h5ad/NOPE/genes")) == {"success": -1, "message": "No dataset found with that ID"}


# ---------------------------------------------------------------- h5ad obs info


def obs_adata():
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["IHC", "OHC", "IHC", "OHC"]),
            "age": pd.Categorical(["P1", None, "P7", "P1"]),
            "has_na": pd.Categorical(["x", "NA", None, "x"]),
            "n_genes": [10, 20, 30, 40],
            "replicate": pd.Categorical(["r1", "r2", "r1", "r2"]),
            "time_point_order": [1, 2, 1, 2],
        },
        index=[f"c{i}" for i in range(4)],
    )
    adata = anndata.AnnData(X=np.ones((4, 2), dtype=np.float32), obs=obs)
    adata.obsm["X_umap"] = np.zeros((4, 2))
    adata.obsm["X_pca"] = np.zeros((4, 2))
    return adata


@pytest.fixture
def h5ad_client(monkeypatch):
    datasets = {"DS1": SINGLE_CELL, "SP1": SPATIAL}
    monkeypatch.setattr(h5ad.geardb, "get_dataset_by_id", lambda dataset_id: datasets.get(dataset_id))
    return make_client(h5ad.H5ad, "/h5ad/<dataset_id>")


def test_h5ad_obs_columns_and_levels(h5ad_client, monkeypatch):
    calls = record_loaders(monkeypatch, h5ad, obs_adata())
    result = ok_json(h5ad_client.get("/h5ad/DS1?analysis_id=A1"))
    assert calls == [("get_adata_shadow", "A1", "DS1", None)]
    assert result["success"] == 1
    assert result["num_obs"] == 4
    # replicate and time_point_order are hidden; embedding coordinates are offered as columns
    assert result["obs_columns"] == ["cell_type", "age", "has_na", "n_genes", "X_umap_1", "X_umap_2", "X_pca_1", "X_pca_2"]
    assert result["has_replicates"] == 1
    # only categorical columns have levels; NA is offered once for columns with missing values
    assert result["obs_levels"] == {"cell_type": ["IHC", "OHC"], "age": ["P1", "P7", "NA"], "has_na": ["NA", "x"]}
    assert result["obs_levels_truncated"] == {}


def test_h5ad_without_replicates_or_embeddings(h5ad_client, monkeypatch):
    obs = pd.DataFrame({"cluster": pd.Categorical(["a", "b"])}, index=["c0", "c1"])
    record_loaders(monkeypatch, h5ad, anndata.AnnData(X=np.ones((2, 1), dtype=np.float32), obs=obs))
    result = ok_json(h5ad_client.get("/h5ad/SP1"))
    assert result["obs_columns"] == ["cluster"]
    assert result["has_replicates"] == 0


def test_h5ad_spatial_uses_spatial_loader(h5ad_client, monkeypatch):
    calls = record_loaders(monkeypatch, h5ad, obs_adata())
    ok_json(h5ad_client.get("/h5ad/SP1"))
    assert [c[0] for c in calls] == ["get_spatial_adata"]


@pytest.mark.parametrize("error, message", [
    (FileNotFoundError("gone"), "No dataset file found."),
    (ValueError("corrupt"), "Encountered an issue with loading the AnnData object"),
])
def test_h5ad_loader_errors(h5ad_client, monkeypatch, error, message):
    record_loaders(monkeypatch, h5ad, error=error)
    assert ok_json(h5ad_client.get("/h5ad/DS1")) == {"success": -1, "message": message}


def test_h5ad_unknown_dataset(h5ad_client):
    assert ok_json(h5ad_client.get("/h5ad/NOPE")) == {"success": -1, "message": "No dataset found with that ID"}


# ---------------------------------------------------------------- HiGlass gene lookup


@pytest.fixture
def higlass_client(monkeypatch):
    """A test client plus the parsed query of each request sent to HiGlass; HiGlass returns client.hits."""
    sent = []

    def fake_get(url, params=None, **kwargs):
        prepared = requests.Request("GET", url, params=params).prepare()
        parts = urlsplit(prepared.url)
        sent.append((f"{parts.scheme}://{parts.netloc}{parts.path}", parse_qs(parts.query)))
        if isinstance(client.hits, Exception):
            raise client.hits
        return types.SimpleNamespace(raise_for_status=lambda: None, json=lambda: client.hits)

    monkeypatch.setattr(higlass.requests, "get", fake_get)
    client = make_client(higlass.HiGlassGene, "/higlass/genes/<gene_symbol>")
    client.sent = sent
    client.hits = []
    return client


def test_higlass_prefers_exact_case_insensitive_match(higlass_client):
    higlass_client.hits = [{"geneName": "Sox21", "chr": "chr14"}, {"geneName": "SOX2", "chr": "chr3", "txStart": 5}]
    result = ok_json(higlass_client.get("/higlass/genes/sox2?assembly=mm10"))
    assert result == {"geneName": "SOX2", "chr": "chr3", "txStart": 5}
    assert higlass_client.sent == [(f"{higlass.HIGLASS_URL}/suggest/",
                                    {"d": [higlass.ANNOTATION_BEDDB_UID["mm10"]], "ac": ["sox2"]})]


def test_higlass_falls_back_to_first_hit(higlass_client):
    higlass_client.hits = [{"geneName": "Sox21"}, {"geneName": "Sox2ot"}]
    assert ok_json(higlass_client.get("/higlass/genes/sox2?assembly=hg38")) == {"geneName": "Sox21"}


def test_higlass_no_hits_or_request_failure_returns_null(higlass_client):
    assert ok_json(higlass_client.get("/higlass/genes/nope?assembly=hg38")) is None
    higlass_client.hits = requests.ConnectionError("down")
    assert ok_json(higlass_client.get("/higlass/genes/sox2?assembly=hg38")) is None


def test_higlass_gene_name_is_url_encoded(higlass_client):
    """A symbol containing & or # stays one ac value instead of adding query parameters."""
    ok_json(higlass_client.get("/higlass/genes/a%26b%23c?assembly=hg38"))
    ((_url, query),) = higlass_client.sent
    assert query == {"d": [higlass.ANNOTATION_BEDDB_UID["hg38"]], "ac": ["a&b#c"]}


@pytest.mark.parametrize("query", ["", "?assembly=galGal6"])
def test_higlass_unsupported_assembly_is_400(higlass_client, query):
    resp = higlass_client.get(f"/higlass/genes/sox2{query}")
    assert resp.status_code == 400
    result = resp.get_json()
    assert result["success"] == -1
    assert "not supported" in result["message"]
    assert higlass_client.sent == []


# ---------------------------------------------------------------- SVG expression data


def write_svg_h5ad(path, with_cell_type=True, symbols=("G1", "G2")):
    obs = pd.DataFrame(index=["t1", "t2", "t3"])
    if with_cell_type:
        obs["cell_type"] = pd.Categorical(["a", "a", "b"])
    var = pd.DataFrame({"gene_symbol": list(symbols)}, index=["E1", "E2"])
    X = np.array([[1, 10], [3, 20], [5, 30]], dtype=np.float32)
    anndata.AnnData(X=X, obs=obs, var=var).write_h5ad(path)


@pytest.fixture
def svg_client(tmp_path, monkeypatch):
    h5_path = tmp_path / "DS1.h5ad"
    dataset = types.SimpleNamespace(id="DS1", get_file_path=lambda: str(h5_path))
    monkeypatch.setattr(svg_data.geardb, "get_dataset_by_id", lambda dataset_id: dataset if dataset_id == "DS1" else None)
    client = make_client(svg_data.SvgData, "/plot/<dataset_id>/svg")
    client.h5_path = h5_path
    return client


def test_svg_scores_and_cell_type_means(svg_client):
    write_svg_h5ad(svg_client.h5_path)
    result = ok_json(svg_client.get("/plot/DS1/svg?gene=G1"))
    assert result["success"] == 1
    assert result["message"] == ""
    assert result["scores"] == {
        "dataset": {"dataset_id": "DS1", "min": 1.0, "max": 30.0},
        "gene": {"gene": "G1", "min": 1.0, "max": 5.0},
        "tissue": {"t1": {"min": 1.0, "max": 10.0}, "t2": {"min": 3.0, "max": 20.0}, "t3": {"min": 5.0, "max": 30.0}},
    }
    assert result["data"] == {"t1": 1.0, "t2": 3.0, "t3": 5.0, "a--mean": 2.0, "b--mean": 5.0}


def test_svg_without_cell_type_returns_tissue_values(svg_client):
    write_svg_h5ad(svg_client.h5_path, with_cell_type=False)
    result = ok_json(svg_client.get("/plot/DS1/svg?gene=G2"))
    assert result["data"] == {"t1": 10.0, "t2": 20.0, "t3": 30.0}


def test_svg_clip_and_user_defined_range(svg_client):
    write_svg_h5ad(svg_client.h5_path)
    result = ok_json(svg_client.get("/plot/DS1/svg?gene=G1&expression_min_clip=2&vmin=0&vmax=4"))
    assert result["scores"]["gene"] == {"gene": "G1", "min": 2.0, "max": 5.0}
    assert result["scores"]["user_defined"] == {"min": 0.0, "max": 4.0}
    # vmax "0" (the client default) means the gene's max
    result = ok_json(svg_client.get("/plot/DS1/svg?gene=G1&vmax=0"))
    assert result["scores"]["user_defined"] == {"min": 0.0, "max": 5.0}


@pytest.mark.parametrize("query", ["vmax=0.0", "vmax=", "vmin=&vmax=0", "vmin=0&vmax=0"])
def test_svg_zero_or_blank_vmax_uses_gene_max(svg_client, query):
    write_svg_h5ad(svg_client.h5_path)
    result = ok_json(svg_client.get(f"/plot/DS1/svg?gene=G1&{query}"))
    assert result["success"] == 1, result["message"]
    if "user_defined" in result["scores"]:
        assert result["scores"]["user_defined"] == {"min": 0.0, "max": 5.0}


@pytest.mark.parametrize("query", ["vmin=abc", "vmax=abc", "vmin=1&vmax=x"])
def test_svg_non_numeric_range_gives_clear_error(svg_client, query):
    write_svg_h5ad(svg_client.h5_path)
    result = ok_json(svg_client.get(f"/plot/DS1/svg?gene=G1&{query}"))
    assert result == {"success": -1, "message": "The vmin and vmax parameters must be numbers."}


def test_svg_max_below_min_gives_clear_error(svg_client):
    write_svg_h5ad(svg_client.h5_path)
    result = ok_json(svg_client.get("/plot/DS1/svg?gene=G1&vmin=4&vmax=2"))
    assert result == {"success": -1, "message": "Maximum cap value must be greater than minimum cap value."}


def test_svg_duplicate_symbol_uses_first_ensembl_id(svg_client):
    write_svg_h5ad(svg_client.h5_path, with_cell_type=False, symbols=("G1", "G1"))
    result = ok_json(svg_client.get("/plot/DS1/svg?gene=G1"))
    assert result["success"] == 2
    assert "Multiple Ensemble IDs" in result["message"]
    assert result["data"] == {"t1": 1.0, "t2": 3.0, "t3": 5.0}


def test_svg_projection_error_is_reported(svg_client, monkeypatch):
    write_svg_h5ad(svg_client.h5_path)

    def fail(adata, dataset_id, projection_id):
        raise PlotError("projection P1 not found")

    monkeypatch.setattr(svg_data, "create_projection_adata", fail)
    result = ok_json(svg_client.get("/plot/DS1/svg?gene=G1&projection_id=P1"))
    assert result == {"success": -1, "message": "projection P1 not found"}


@pytest.mark.parametrize("url, message", [
    ("/plot/DS1/svg", "Request needs both dataset id and gene symbol."),
    ("/plot/DS1/svg?gene=G1&expression_min_clip=abc", "The expression_min_clip parameter must be a number."),
    ("/plot/NOPE/svg?gene=G1", "The dataset was not found."),
    ("/plot/DS1/svg?gene=NOPE", "The gene symbol 'NOPE' was not found in the dataset."),
])
def test_svg_error_messages(svg_client, url, message):
    write_svg_h5ad(svg_client.h5_path)
    assert ok_json(svg_client.get(url)) == {"success": -1, "message": message}


def test_svg_missing_file(svg_client):
    assert ok_json(svg_client.get("/plot/DS1/svg?gene=G1")) == {"success": -1, "message": "The h5ad file was not found."}


def test_svg_missing_gene_symbol_column(svg_client):
    anndata.AnnData(X=np.ones((2, 2), dtype=np.float32)).write_h5ad(svg_client.h5_path)
    result = ok_json(svg_client.get("/plot/DS1/svg?gene=G1"))
    assert result == {"success": -1, "message": "The h5ad is missing the gene_symbol column."}


# ---------------------------------------------------------------- dataset displays


@pytest.fixture
def display_client(monkeypatch):
    display = {"id": 3, "dataset_id": "DS1", "user_id": 5, "label": "UMAP", "plot_type": "umap_static",
               "plotly_config": {"colorize_legend_by": "cell_type"}}
    monkeypatch.setattr(dataset_display.geardb, "get_display_by_id",
                        lambda display_id: display if str(display_id) == "3" else {})
    return make_client(dataset_display.DatasetDisplay, "/displays/<display_id>")


def test_dataset_display_found(display_client):
    result = ok_json(display_client.get("/displays/3"))
    assert result["label"] == "UMAP"
    assert result["plotly_config"] == {"colorize_legend_by": "cell_type"}


def test_dataset_display_not_found_is_404(display_client):
    resp = display_client.get("/displays/99")
    assert resp.status_code == 404
    assert resp.get_json() == {"message": "Display not found"}
