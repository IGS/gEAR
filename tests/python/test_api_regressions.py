"""
Regression tests for closed GitHub issues in the Flask API resources (www/api/resources).

Each test names its issue. Resources are called through Flask's test client with small
in-memory AnnData objects standing in for the dataset loaders.
"""

import base64
import importlib
import sys
import types

import anndata
import h5py
import numpy as np
import pandas as pd
import pytest
from flask import Flask
from flask_restful import Api
from resources import available_display_types, common, h5ad, plotly_data

SINGLE_CELL = types.SimpleNamespace(id="DS1", dtype="single-cell-rnaseq")


def make_client(resource, route):
    app = Flask(__name__)
    Api(app).add_resource(resource, route)
    return app.test_client()


def ok_json(resp):
    assert resp.status_code == 200, resp.get_data(as_text=True)
    return resp.get_json()


# ---------------------------------------------------------------- multigene plots (mg_plotly_data)

@pytest.fixture(scope="module")
def mg_plotly_data():
    """Import resources.mg_plotly_data, standing in for diffxpy (only used for differential expression) if it is missing."""
    with pytest.MonkeyPatch.context() as mp:
        try:
            import diffxpy.api  # noqa: F401
        except ImportError:
            diffxpy = types.ModuleType("diffxpy")
            diffxpy.api = types.ModuleType("diffxpy.api")
            mp.setitem(sys.modules, "diffxpy", diffxpy)
            mp.setitem(sys.modules, "diffxpy.api", diffxpy.api)
        mp.delitem(sys.modules, "gear.mg_plotting", raising=False)
        mp.delitem(sys.modules, "resources.mg_plotly_data", raising=False)
        yield importlib.import_module("resources.mg_plotly_data")


def build_mg_adata():
    """4 cells: two in cluster "a" and two in cluster "b", each with a different condition."""
    obs = pd.DataFrame(
        {
            "cluster": pd.Categorical(["a", "a", "b", "b"]),
            "condition": pd.Categorical(["x", "y", "x", "y"]),
        },
        index=["c0", "c1", "c2", "c3"],
    )
    var = pd.DataFrame({"gene_symbol": ["G1", "G2"]}, index=pd.Index(["E1", "E2"], name="ensembl_id"))
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]], dtype=np.float32)
    return anndata.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def mg_client(mg_plotly_data, monkeypatch):
    mod = mg_plotly_data
    monkeypatch.setattr(mod.geardb, "get_dataset_by_id", lambda *a, **k: SINGLE_CELL)
    monkeypatch.setattr(mod, "get_adata_from_analysis", lambda *a, **k: build_mg_adata())
    return make_client(mod.MGPlotlyData, "/plot/<dataset_id>/mg_plotly")


def mg_plot(client, **body):
    body = {"gene_symbols": ["G1", "G2"], "plot_type": "heatmap", "primary_col": "cluster", **body}
    return ok_json(client.post("/plot/DS1/mg_plotly", json=body))


def decode(values):
    """A trace's x/y/z values as a list or array, decoding Plotly's base64 typed arrays."""
    if isinstance(values, dict):
        shape = [int(n) for n in str(values.get("shape", "-1")).split(",")]
        return np.frombuffer(base64.b64decode(values["bdata"]), dtype=values["dtype"]).reshape(shape)
    return values


def heatmap_z(result):
    """The main heatmap's z matrix (genes x observation columns)."""
    (trace,) = [t for t in result["plot_json"]["data"] if t["type"] == "heatmap"]
    return np.asarray(decode(trace["z"]), dtype=float)


def test_heatmap_without_matrixplot_keeps_one_column_per_observation(mg_client):
    """Regression for #1324: with matrixplot off, cells sharing a group are separate heatmap columns, not averaged."""
    result = mg_plot(mg_client, matrixplot=False)
    assert result["success"] == 1, result["message"]
    z = heatmap_z(result)
    assert z.shape == (2, 4)
    # log2(expression + 1) of each of the four cells for G1
    assert sorted(z[0]) == pytest.approx(np.log2([2, 4, 6, 8]))


def test_heatmap_matrixplot_averages_observations_within_a_group(mg_client):
    """Regression for #1324: with matrixplot on, cells in the same group are averaged into one column."""
    result = mg_plot(mg_client, matrixplot=True)
    assert result["success"] == 1, result["message"]
    z = heatmap_z(result)
    assert z.shape == (2, 2)
    # G1 means: a = (1 + 3) / 2, b = (5 + 7) / 2
    assert z[0] == pytest.approx(np.log2([2 + 1, 6 + 1]))


@pytest.mark.parametrize(
    "body, expected",
    [
        ({"obs_filters": {"no_such_column": ["a"]}}, "Filter series no_such_column not found"),
        ({"primary_col": "no_such_column", "sort_order": {"no_such_column": ["a"]}},
         "Sort order series no_such_column not found"),
    ],
)
def test_missing_filter_or_sort_column_gives_clear_error(mg_client, body, expected):
    """Regression for #919: a filter or sort column missing from obs returns a clear error, not a 500."""
    result = mg_plot(mg_client, **body)
    assert result["success"] == -1
    assert expected in result["message"]


def test_stale_sort_level_gives_clear_error(mg_client):
    """Regression for #919: a saved order naming a level the dataset no longer has gives a clear error, not pandas' "new_categories" error."""
    result = mg_plot(
        mg_client, sort_order={"cluster": ["b", "a", "renamed"]}, obs_filters={"cluster": ["a", "b"]}
    )
    assert result["success"] == -1
    assert "At least one value for filter series cluster" in result["message"]
    assert "new_categories" not in result["message"]


@pytest.mark.parametrize("plot_type", ["heatmap", "dotplot"])
def test_sorted_group_without_a_filter_on_it_still_plots(mg_client, plot_type):
    """Regression for #919: sorting the primary group while filtering on a different column plots instead of raising KeyError (500)."""
    result = mg_plot(
        mg_client, plot_type=plot_type, sort_order={"cluster": ["b", "a"]}, obs_filters={"condition": ["x"]}
    )
    assert result["success"] == 1, result["message"]
    if plot_type == "heatmap":
        # One cell per cluster is left (c0 and c2), plotted in the saved order b, a
        assert heatmap_z(result)[0] == pytest.approx(np.log2([5 + 1, 1 + 1]))


@pytest.mark.parametrize("categorical", [True, False])
def test_filter_on_na_level_selects_cells_with_missing_values(mg_client, mg_plotly_data, monkeypatch, categorical):
    """The curator offers "NA" as a level for columns with missing values, so filtering on it must select those cells, not raise (500)."""
    adata = build_mg_adata()
    condition = ["x", None, "x", None]
    adata.obs["condition"] = pd.Categorical(condition) if categorical else pd.Series(condition, index=adata.obs.index, dtype=object)
    monkeypatch.setattr(mg_plotly_data, "get_adata_from_analysis", lambda *a, **k: adata)

    result = mg_plot(mg_client, matrixplot=False, obs_filters={"condition": ["NA"]})
    assert result["success"] == 1, result["message"]
    # Only c1 and c3 (the cells with no condition) are left
    assert sorted(heatmap_z(result)[0]) == pytest.approx(np.log2([3 + 1, 7 + 1]))

    result = mg_plot(mg_client, matrixplot=False, obs_filters={"condition": ["x", "NA"]})
    assert result["success"] == 1, result["message"]
    assert heatmap_z(result).shape == (2, 4)


@pytest.mark.xfail(strict=True, reason=(
    "mg_plotly_data checks the saved order against the categories left after all filters, so a level "
    "emptied by another filter is reported as missing from the dataset ('Please update curation')"
))
def test_level_emptied_by_another_filter_is_not_reported_as_stale(mg_client, mg_plotly_data, monkeypatch):
    """Regression for #919: only levels missing from the dataset are a stale curation; a level another filter empties is not."""
    adata = build_mg_adata()
    adata.obs["condition"] = pd.Categorical(["x", "x", "y", "y"])     # cluster "b" has no "x" cells
    monkeypatch.setattr(mg_plotly_data, "get_adata_from_analysis", lambda *a, **k: adata)
    result = mg_plot(
        mg_client, sort_order={"cluster": ["b", "a"]}, obs_filters={"cluster": ["b", "a"], "condition": ["x"]}
    )
    assert result["success"] == 1, result["message"]


# ---------------------------------------------------------------- single-gene plots (plotly_data)

@pytest.mark.parametrize("name", ["PC1", "condition"])
def test_gene_named_like_an_obs_column_still_plots(monkeypatch, name):
    """Regression for #878: a gene (e.g. a projection pattern) named like an obs column plots instead of clashing with that column."""
    obs = pd.DataFrame(
        {name: pd.Categorical(["p", "p", "q", "q"]) if name == "condition" else [0.5, 1.5, 2.5, 3.5]},
        index=[f"cell{i}" for i in range(4)],
    )
    var = pd.DataFrame({"gene_symbol": [name, "Other"]}, index=[name, "Other"])
    X = np.array([[10, 0], [20, 0], [30, 0], [40, 0]], dtype=np.float32)
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    monkeypatch.setattr(plotly_data.geardb, "get_dataset_by_id", lambda *a, **k: SINGLE_CELL)
    monkeypatch.setattr(plotly_data, "get_adata_from_analysis", lambda *a, **k: adata)
    client = make_client(plotly_data.PlotlyData, "/plot/<dataset_id>")

    result = ok_json(client.post("/plot/DS1", json={
        "gene_symbol": name, "plot_type": "scatter", "x_axis": name, "y_axis": "raw_value",
    }))
    assert result["success"] == 1, result["message"]
    points = sorted(
        (x, float(y))
        for trace in result["plot_json"]["data"]
        for x, y in zip(decode(trace["x"]), decode(trace["y"]))
    )
    expected_x = obs[name].tolist()
    assert points == sorted(zip(expected_x, [10.0, 20.0, 30.0, 40.0]))


# ---------------------------------------------------------------- available display types

def display_types(monkeypatch, tmp_path, adata, resource="AvailableDisplayTypes"):
    ds = types.SimpleNamespace(
        id="DS1", dtype="single-cell-rnaseq", get_file_path=lambda: str(tmp_path / "datasets" / "DS1.h5ad")
    )
    mod = available_display_types
    monkeypatch.setattr(mod.geardb, "get_dataset_by_id", lambda *a, **k: ds)
    monkeypatch.setattr(mod, "get_adata_shadow", lambda *a, **k: adata)
    client = make_client(getattr(mod, resource), "/h5ad/<dataset_id>/types")
    return ok_json(client.post("/h5ad/DS1/types", json={"dataset_id": "DS1"}))


def embedding_adata(obs_columns, obsm_key=None):
    n = 4
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(n)])
    obs["cluster"] = pd.Categorical(["a", "a", "b", "b"])
    for col, kind in obs_columns.items():
        obs[col] = np.linspace(0, 1, n) if kind == "numeric" else pd.Categorical(["u", "v", "u", "v"])
    adata = anndata.AnnData(X=np.ones((n, 2), dtype=np.float32), obs=obs)
    if obsm_key:
        adata.obsm[obsm_key] = np.zeros((n, 2))
    return adata


def test_static_embeddings_from_numeric_columns_with_any_name(monkeypatch, tmp_path):
    """Regression for #876: two numeric obs columns enable static tSNE/UMAP/PCA plots even if their names don't mention tsne or umap."""
    adata = embedding_adata({"Embeddings_1": "numeric", "Embeddings_2": "numeric"})
    types_ = display_types(monkeypatch, tmp_path, adata)
    assert types_["tsne_static"] and types_["umap_static"] and types_["pca_static"]
    assert types_["tsne/umap_dynamic"]


@pytest.mark.parametrize(
    "obsm_key, static_type", [("X_tsne", "tsne_static"), ("X_umap", "umap_static"), ("X_pca", "pca_static")]
)
def test_static_embedding_from_obsm(monkeypatch, tmp_path, obsm_key, static_type):
    """Regression for #876: an embedding in obsm enables its static plot, without any obs column naming it."""
    types_ = display_types(monkeypatch, tmp_path, embedding_adata({}, obsm_key))
    assert types_[static_type] is True
    assert types_["tsne/umap_dynamic"] is True


def test_embedding_like_column_name_alone_does_not_enable_static_plots(monkeypatch, tmp_path):
    """Regression for #876: availability no longer depends on column names, so a categorical "tSNE_group" enables nothing."""
    types_ = display_types(monkeypatch, tmp_path, embedding_adata({"tSNE_group": "categorical", "UMAP_label": "categorical"}))
    assert not types_["tsne_static"] and not types_["umap_static"] and not types_["pca_static"]
    assert not types_["tsne/umap_dynamic"]


# ---------------------------------------------------------------- obs levels (h5ad)

def test_high_cardinality_categorical_column_reported_as_truncated(monkeypatch):
    """Regression for #1320: a categorical column with more than 50 levels is listed under obs_levels_truncated, not dropped."""
    n = 60
    obs = pd.DataFrame(
        {
            "barcode": pd.Categorical([f"BC{i:02d}" for i in range(n)]),
            "fifty": pd.Categorical([f"L{i % 50:02d}" for i in range(n)]),
            "cluster": pd.Categorical(["a", "b", "c"] * (n // 3)),
        },
        index=[f"cell{i}" for i in range(n)],
    )
    adata = anndata.AnnData(X=np.ones((n, 1), dtype=np.float32), obs=obs)
    monkeypatch.setattr(h5ad.geardb, "get_dataset_by_id", lambda *a, **k: SINGLE_CELL)
    monkeypatch.setattr(h5ad, "get_adata_shadow", lambda *a, **k: adata)

    result = ok_json(make_client(h5ad.H5ad, "/h5ad/<dataset_id>").get("/h5ad/DS1"))
    assert result["success"] == 1
    assert "barcode" in result["obs_columns"]
    assert result["obs_levels_truncated"] == {"barcode": [f"BC{i:02d}" for i in range(n)]}
    assert set(result["obs_levels"]) == {"fifty", "cluster"}     # exactly 50 levels is not truncated
    assert result["obs_levels"]["cluster"] == ["a", "b", "c"]


# ---------------------------------------------------------------- projections (common.create_projection_adata)

def test_projection_built_twice_from_backed_file_releases_it(monkeypatch, tmp_path):
    """Regression for #1062 / #927: building a projection closes the backed dataset file, so it can be built again (no BlockingIOError)."""
    cells = [f"cell{i}" for i in range(4)]
    obs = pd.DataFrame({"cluster": pd.Categorical(["a", "a", "b", "b"])}, index=cells)
    dataset_path = tmp_path / "DS1.h5ad"
    anndata.AnnData(X=np.ones((4, 3), dtype=np.float32), obs=obs).write_h5ad(dataset_path)

    projection_dir = tmp_path / "projections" / "by_dataset" / "DS1"
    projection_dir.mkdir(parents=True)
    patterns = pd.DataFrame({"P1": [0.1, 0.2, 0.3, 0.4], "P2": [1.0, 2.0, 3.0, 4.0]}, index=cells)
    patterns.to_csv(projection_dir / "proj1.csv")
    monkeypatch.setattr(common, "PROJECTIONS_BASE_DIR", tmp_path / "projections")

    for _ in range(2):
        dataset_adata = anndata.read_h5ad(dataset_path, backed="r")
        projection = common.create_projection_adata(dataset_adata, "DS1", "proj1")

        assert not dataset_adata.file.is_open
        assert projection.var["gene_symbol"].tolist() == ["P1", "P2"]
        assert np.allclose(projection.X, patterns.to_numpy())
        assert projection.obs["cluster"].tolist() == ["a", "a", "b", "b"]
        # The projection is in memory: it holds no file that could block the next request
        assert not projection.isbacked

    # Nothing still holds the dataset file: opening it for writing would fail otherwise
    h5py.File(dataset_path, "r+").close()
    # No projection h5ad is left behind next to the CSV
    assert sorted(p.name for p in projection_dir.iterdir()) == ["proj1.csv"]
