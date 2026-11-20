"""
Helpers of the spatial plot resources: spatial_scanpy_data (core) and spatialpanel, plus the
SpatialPanel endpoint's cache paths. spatialpanel imports SpatialData, so its tests are marked spatial.
"""

import importlib
import types

import anndata
import numpy as np
import pandas as pd
import pytest
from flask import Flask
from flask_restful import Api
from resources import spatial_scanpy_data


@pytest.fixture(scope="module")
def spatialpanel():
    pytest.importorskip("spatialdata")
    pytest.importorskip("bokeh")
    return importlib.import_module("resources.spatialpanel")


def build_adata(clusters, gene_symbols=("Actb", "Gfap")):
    n_cells = len(clusters)
    rng = np.random.default_rng(0)
    obs = pd.DataFrame(
        {
            "clusters": pd.Categorical(clusters),
            "spatial1": np.arange(n_cells, dtype=float),
            "spatial2": np.arange(n_cells, dtype=float) * 2,
            "n_genes_by_counts": np.full(n_cells, 100),
        },
        index=[f"cell{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(
        {"gene_symbol": list(gene_symbols)}, index=[f"ENSG{i}" for i in range(len(gene_symbols))]
    )
    X = rng.uniform(size=(n_cells, len(gene_symbols))).astype(np.float32)
    return anndata.AnnData(X=X, obs=obs, var=var)


def cluster_colors(adata):
    """Colors per cluster as scanpy reads them: uns colors paired with the categories by position."""
    return dict(zip(adata.obs["clusters"].cat.categories, adata.uns["clusters_colors"]))


# ---------------------------------------------------------------- spatial_scanpy_data


def test_scanpy_map_colors_sorts_numeric_clusters():
    adata = spatial_scanpy_data.map_colors(build_adata(["1", "10", "2", "1"]))
    alphabet = spatial_scanpy_data.px.colors.qualitative.Alphabet
    assert adata.obs["clusters"].cat.categories.tolist() == ["1", "2", "10"]
    assert cluster_colors(adata) == {"1": alphabet[0], "2": alphabet[1], "10": alphabet[2]}


def test_scanpy_map_colors_accepts_named_clusters():
    adata = spatial_scanpy_data.map_colors(build_adata(["T cell", "B cell", "T cell"]))
    alphabet = spatial_scanpy_data.px.colors.qualitative.Alphabet
    assert cluster_colors(adata) == {"B cell": alphabet[0], "T cell": alphabet[1]}


def test_scanpy_map_colors_keeps_each_clusters_own_color():
    adata = build_adata(["1", "10", "2", "10"])
    adata.obs["colors"] = ["#111111", "#aaaaaa", "#222222", "#aaaaaa"]
    adata = spatial_scanpy_data.map_colors(adata)
    assert cluster_colors(adata) == {"1": "#111111", "2": "#222222", "10": "#aaaaaa"}


def test_normalize_gene_symbols_is_case_insensitive_and_drops_unknown_genes():
    dataset_genes = ["Actb", "Gfap", "Sox2"]
    assert spatial_scanpy_data.normalize_gene_symbols(dataset_genes, ["SOX2", "nope", "actb"]) == ["Sox2", "Actb"]
    assert spatial_scanpy_data.normalize_gene_symbols(dataset_genes, []) == []


@pytest.mark.xfail(
    strict=True,
    raises=(AttributeError, NameError),
    reason="create_plot crashes under scanpy 1.12: one gene calls .remove() on the violin legend scanpy no "
    "longer draws; several genes use an undefined name 'clusters'. No page calls /spatial_scanpy.",
)
@pytest.mark.parametrize("genes", [["Actb"], ["Actb", "Gfap"]])
def test_scanpy_create_plot(genes):
    adata = build_adata(["1", "2", "3"] * 4, gene_symbols=genes)
    adata.var_names = genes
    rng = np.random.default_rng(1)
    adata.obsm["spatial"] = rng.normal(size=(adata.n_obs, 2))
    adata.obsm["X_umap"] = rng.normal(size=(adata.n_obs, 2))
    adata.uns["has_images"] = False
    spatial_scanpy_data.create_plot(spatial_scanpy_data.map_colors(adata), genes, False)


@pytest.fixture
def scanpy_client(monkeypatch):
    app = Flask(__name__)
    Api(app).add_resource(spatial_scanpy_data.SpatialScanpyData, "/plot/<dataset_id>/spatial_scanpy")
    return app.test_client()


def test_scanpy_endpoint_requires_genes(scanpy_client):
    data = scanpy_client.post("/plot/DS1/spatial_scanpy", json={"gene_symbols": []}).get_json()
    assert data == {"success": -1, "message": "Request needs at least one gene symbol."}


def test_scanpy_endpoint_reports_loader_failure(scanpy_client, monkeypatch):
    def fail(*args, **kwargs):
        raise FileNotFoundError("gone")

    monkeypatch.setattr(spatial_scanpy_data, "get_spatial_adata", fail)
    data = scanpy_client.post("/plot/DS1/spatial_scanpy", json={"gene_symbols": ["Actb"]}).get_json()
    assert data == {"success": -1, "message": "Could not create dataset object using analysis."}


def test_scanpy_endpoint_reports_projection_errors(scanpy_client, monkeypatch):
    def fail(adata, dataset_id, projection_id):
        raise spatial_scanpy_data.PlotError("Projection P1 is not available")

    monkeypatch.setattr(spatial_scanpy_data, "get_spatial_adata", lambda *a, **k: build_adata(["1"]))
    monkeypatch.setattr(spatial_scanpy_data, "create_projection_adata", fail)
    data = scanpy_client.post(
        "/plot/DS1/spatial_scanpy", json={"gene_symbols": ["Actb"], "projection_id": "P1"}
    ).get_json()
    assert data == {"success": -1, "message": "Projection P1 is not available"}


# ---------------------------------------------------------------- spatialpanel helpers


@pytest.mark.spatial
@pytest.mark.parametrize(
    "clusters, expected",
    [(["10", "2", "1"], ["1", "2", "10"]), ([3, 1, 2], [1, 2, 3]), (["b", "10", "a"], ["10", "a", "b"])],
)
def test_sort_clusters(spatialpanel, clusters, expected):
    assert spatialpanel.sort_clusters(clusters) == expected


@pytest.mark.spatial
def test_normalize_searched_gene(spatialpanel):
    genes = {"Actb", "Gfap"}
    assert spatialpanel.normalize_searched_gene(genes, "ACTB") == "Actb"
    assert spatialpanel.normalize_searched_gene(genes, "gfap") == "Gfap"
    assert spatialpanel.normalize_searched_gene(genes, "Sox2") is None


@pytest.mark.spatial
def test_normalize_image_array_stretches_to_uint8_rgb(spatialpanel):
    arr = np.linspace(0, 1000, 10_000, dtype=np.float64).reshape(100, 100)
    arr[0, 0] = 1e9  # outlier, clipped at the 99th percentile
    img = spatialpanel.normalize_image_array(arr)
    assert img.shape == (100, 100, 3)
    assert img.dtype == np.uint8
    assert img.min() == 0
    assert img.max() >= 254
    assert np.array_equal(img[..., 0], img[..., 2])


@pytest.mark.spatial
def test_normalize_image_array_keeps_uint8_values(spatialpanel):
    gray = np.arange(12, dtype=np.uint8).reshape(3, 4, 1)
    assert np.array_equal(spatialpanel.normalize_image_array(gray), np.repeat(gray, 3, axis=-1))
    rgb = np.zeros((2, 2, 3), dtype=np.uint8)
    assert spatialpanel.normalize_image_array(rgb) is rgb


@pytest.mark.spatial
def test_normalize_image_array_rejects_missing_image(spatialpanel):
    with pytest.raises(ValueError, match="Image array is None"):
        spatialpanel.normalize_image_array(None)


@pytest.mark.spatial
@pytest.mark.parametrize("found_img, is_cool, palette", [(False, False, "glasbey_light"),
                                                         (True, False, "b_glasbey_bw"),
                                                         (True, True, "glasbey_cool")])
def test_panel_map_colors_assigns_palette_in_cluster_order(spatialpanel, found_img, is_cool, palette):
    df = pd.DataFrame({"clusters": pd.Categorical(["10", "2", "1", "2"])})
    df = spatialpanel.map_colors(df, found_img, is_cool)
    swatch = getattr(spatialpanel.cc, palette)
    assert dict(zip(df["clusters"], df["colors"])) == {"1": swatch[0], "2": swatch[1], "10": swatch[2]}


@pytest.mark.spatial
def test_panel_map_colors_keeps_existing_colors(spatialpanel):
    """A cluster category with no rows used to raise IndexError when a colors column was present."""
    df = pd.DataFrame({
        "clusters": pd.Categorical(["1", "2"], categories=["1", "2", "3"]),
        "colors": ["#111111", "#222222"],
    })
    df = spatialpanel.map_colors(df, False, False)
    assert df["colors"].tolist() == ["#111111", "#222222"]


@pytest.mark.spatial
def test_create_gene_df(spatialpanel):
    adata = build_adata(["2", "1", None, "2"])
    adata.obsm["X_umap"] = np.array([[0.0, 1.0], [2.0, 3.0], [4.0, 5.0], [6.0, 7.0]])
    df = spatialpanel.create_gene_df(adata, "GFAP")
    assert list(df.columns) == ["raw_value", "spatial1", "spatial2", "n_genes_by_counts",
                                "UMAP1", "UMAP2", "clusters", "clusters_cat_codes"]
    assert df.index.tolist() == ["cell0", "cell1", "cell3"]  # NA cluster dropped
    assert np.allclose(df["raw_value"], adata.X[[0, 1, 3], 1])
    assert df["UMAP2"].tolist() == [1.0, 3.0, 7.0]
    assert df["clusters"].tolist() == ["2", "1", "2"]


@pytest.mark.spatial
def test_create_gene_df_without_umap(spatialpanel):
    df = spatialpanel.create_gene_df(build_adata(["1", "2"]), "Actb")
    assert "UMAP1" not in df


@pytest.mark.spatial
def test_create_gene_df_errors(spatialpanel):
    with pytest.raises(ValueError, match="Gene 'Sox2' not found"):
        spatialpanel.create_gene_df(build_adata(["1"]), "Sox2")

    adata = build_adata(["1"])
    del adata.obs["clusters"]
    with pytest.raises(ValueError, match="No cluster information"):
        spatialpanel.create_gene_df(adata, "Actb")

    with pytest.raises(ValueError, match="No observations with valid cluster information"):
        spatialpanel.create_gene_df(build_adata([None, None]), "Actb")


# ---------------------------------------------------------------- SpatialPanel endpoint


@pytest.fixture
def panel(spatialpanel, monkeypatch, tmp_path):
    """SpatialPanel with its cache under tmp_path/cache and a dataset without images."""
    state = types.SimpleNamespace(cache=tmp_path / "cache", handler_calls=[], adata=build_adata(["1", "2", "1"]))
    state.cache.mkdir()

    def initialize_spatial_handler(dataset_id):
        state.handler_calls.append(dataset_id)
        return types.SimpleNamespace(has_images=False)

    monkeypatch.setattr(spatialpanel, "PANEL_CSV_CACHE_DIR", state.cache)
    monkeypatch.setattr(spatialpanel, "initialize_spatial_handler", initialize_spatial_handler)
    monkeypatch.setattr(spatialpanel, "get_table_adata", lambda spatial_obj: state.adata)
    monkeypatch.setattr(spatialpanel, "get_viewer_script", lambda is_zoomed, config: f"<script {config['filename']}>")

    app = Flask(__name__)
    Api(app).add_resource(spatialpanel.SpatialPanel, "/plot/<dataset_id>/spatialpanel")
    client = app.test_client()

    def post(dataset_id="DS1", **body):
        resp = client.post(f"/plot/{dataset_id}/spatialpanel", json=body)
        assert resp.status_code == 200, resp.get_data(as_text=True)
        return resp.get_json()

    state.post = post
    return state


@pytest.mark.spatial
def test_spatial_panel_caches_gene_csv(panel):
    data = panel.post(gene_symbol="actb")
    assert data["success"] == 1, data
    assert data["script"] == "<script actb.csv>"
    df = pd.read_csv(panel.cache / "DS1" / "actb.csv")
    assert df["clusters"].tolist() == [1, 2, 1]
    assert df["colors"].notna().all()

    data = panel.post(gene_symbol="actb")
    assert data["success"] == 1
    assert "Using cached file." in data["message"]


@pytest.mark.spatial
@pytest.mark.parametrize(
    "dataset_id, body, message",
    [
        ("..", {"gene_symbol": "Actb"}, "Invalid dataset ID."),
        ("DS 1", {"gene_symbol": "Actb"}, "Invalid dataset ID."),
        ("DS1", {"gene_symbol": "Actb", "projection_id": "../P1"}, "Invalid projection ID."),
        ("DS1", {"gene_symbol": "../OTHER/Actb"}, "Invalid gene symbol."),
        ("DS1", {}, "Request needs a gene symbol."),
    ],
)
def test_spatial_panel_rejects_unsafe_cache_paths(panel, dataset_id, body, message):
    # Another dataset's cached CSV that "../OTHER/Actb" would otherwise be served from
    (panel.cache / "OTHER").mkdir()
    (panel.cache / "OTHER" / "Actb.csv").write_text("raw_value\n1\n")

    data = panel.post(dataset_id, **body)
    assert data == {"script": None, "success": 0, "message": message}
    assert panel.handler_calls == []
    assert sorted(p.name for p in panel.cache.iterdir()) == ["OTHER"]
