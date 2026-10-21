"""tsne_data helpers and the TSNEData/MGTSNEData endpoints' validation and error responses."""

import types

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import scipy.sparse
from flask import Flask
from flask_restful import Api
from gear.plotting import PlotError
from resources import tsne_data


def build_adata(n_cells=6, gene_symbols=("G1", "G2"), X=None):
    rng = np.random.default_rng(0)
    obs = pd.DataFrame(
        {
            "cluster": pd.Categorical(["a", "b", "c"] * (n_cells // 3), categories=["a", "b", "c"]),
            "tSNE_1": rng.normal(size=n_cells),
            "tSNE_2": rng.normal(size=n_cells),
        },
        index=[f"cell{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(
        {"gene_symbol": list(gene_symbols)}, index=[f"ENSG{i}" for i in range(len(gene_symbols))]
    )
    if X is None:
        X = rng.uniform(size=(n_cells, len(gene_symbols))).astype(np.float32)
    return anndata.AnnData(X=X, obs=obs, var=var)


def analysis(tmp_path, analysis_id="DS1"):
    return types.SimpleNamespace(dataset_id="DS1", id=analysis_id, dataset_path=str(tmp_path / "DS1.h5ad"))


class NonPrimarySpatialAnalysis(tsne_data.SpatialAnalysis):
    """A saved spatial analysis, whose dataset_path is an .h5ad rather than the primary .zarr store."""

    def __init__(self, dataset_path):
        self.id = "A1"
        self.dataset_id = "DS1"
        self._dataset_path = dataset_path

    @property
    def dataset_path(self):
        return self._dataset_path


# ---------------------------------------------------------------- small helpers


@pytest.mark.parametrize(
    "expression_palette, expected",
    [("YlOrRd", "viridis"), ("cividis_r", "viridis"), ("viridis", "plasma"), ("viridis_r", "plasma")],
)
def test_annotation_colormap_differs_from_expression_palette(expression_palette, expected):
    assert tsne_data.annotation_colormap(expression_palette) == expected


@pytest.mark.parametrize("groups, expected", [(1, 1), (16, 1), (17, 2), (32, 2), (33, 3)])
def test_calculate_num_legend_cols(groups, expected):
    assert tsne_data.calculate_num_legend_cols(groups) == expected


@pytest.mark.parametrize(
    "series, expected",
    [
        (pd.Series(pd.Categorical(["a", "b"])), True),
        (pd.Series(["a", "b"]), False),
        (pd.Series([1.0, 2.0]), False),
    ],
)
def test_is_categorical(series, expected):
    assert tsne_data.is_categorical(series) is expected


# ---------------------------------------------------------------- sort_legend


def legend_axes(labels):
    fig, ax = plt.subplots()
    for label in labels:
        ax.scatter([], [], label=label)
    return fig, ax


def test_sort_legend_keeps_category_order_for_order_of_appearance():
    """The endpoint passes obs[colorize_by].unique(); the legend stays in category order (scanpy's)."""
    fig, ax = legend_axes(["a", "b", "c"])
    handles, labels = tsne_data.sort_legend(ax, ["c", "a", "b"], 3)
    plt.close(fig)
    assert labels == ["a", "b", "c"]
    assert [h.get_label() for h in handles] == labels


def test_sort_legend_returns_legend_unchanged_when_sort_order_is_longer():
    fig, ax = legend_axes(["a", "b"])
    handles, labels = tsne_data.sort_legend(ax, ["a", "b", "c"], 3)
    plt.close(fig)
    assert labels == ["a", "b"]


def test_sort_legend_ignores_names_in_sort_order():
    """
    Names in sort_order are intentionally not used: sort_legend keeps the first len(sort_order) legend
    entries in category order. The endpoint passes obs[colorize_by].unique() (order of appearance),
    and category order is what users get until someone asks for something else.
    """
    fig, ax = legend_axes(["a", "b", "c"])
    _, labels = tsne_data.sort_legend(ax, ["x", "y", "z"], 3)     # names not in the legend at all
    _, shorter = tsne_data.sort_legend(ax, ["c", "b"], 3)
    plt.close(fig)
    assert labels == ["a", "b", "c"]
    assert shorter == ["a", "b"]


@pytest.mark.parametrize("n_items, num_cols", [(5, 2), (12, 12), (14, 12), (25, 12), (26, 12)])
def test_horizontal_legend_reads_left_to_right(n_items, num_cols):
    """Matplotlib fills legend columns top to bottom; sort_legend reorders so rows read in order."""
    names = [f"g{i:02d}" for i in range(n_items)]
    fig, ax = legend_axes(names)
    handles, labels = tsne_data.sort_legend(ax, names, num_cols, horizontal_legend=True)
    legend = fig.legend(handles=handles, labels=labels, ncol=num_cols)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    positions = {}
    for text in legend.get_texts():
        box = text.get_window_extent(renderer)
        positions[text.get_text()] = (-round(box.y0), box.x0)
    plt.close(fig)
    assert sorted(positions, key=positions.get) == names


# ---------------------------------------------------------------- validate_args


@pytest.fixture
def validation(monkeypatch):
    """Stub the dataset and analysis lookups behind validate_args; returns a record of the calls."""
    state = types.SimpleNamespace(
        datasets={"DS1": types.SimpleNamespace(id="DS1", dtype="single-cell-rnaseq"),
                  "SP1": types.SimpleNamespace(id="SP1", dtype="spatial")},
        adata=build_adata(),
        analysis_error=None,
        adata_error=None,
        get_analysis_calls=[],
        get_adata_calls=[],
    )

    def get_analysis(analysis, dataset_id, session_id, is_spatial=False):
        state.get_analysis_calls.append({"analysis": analysis, "dataset_id": dataset_id, "is_spatial": is_spatial})
        if state.analysis_error:
            raise state.analysis_error
        ana = types.SimpleNamespace(dataset_id=dataset_id, id=dataset_id, dataset_path="/nonexistent/DS1.h5ad")

        def get_adata(**kwargs):
            state.get_adata_calls.append(kwargs)
            if state.adata_error:
                raise state.adata_error
            return state.adata

        ana.get_adata = get_adata
        return ana

    monkeypatch.setattr(tsne_data.geardb, "get_dataset_by_id", lambda dataset_id: state.datasets.get(dataset_id))
    monkeypatch.setattr(tsne_data, "get_analysis", get_analysis)
    return state


@pytest.mark.parametrize(
    "dataset_id, gene_symbol, message",
    [(None, "G1", "dataset id"), ("DS1", None, "gene symbol"), ("DS1", [], "gene symbol"),
     ("missing", "G1", "No dataset found")],
)
def test_validate_args_rejects_missing_inputs(validation, dataset_id, gene_symbol, message):
    result = tsne_data.validate_args(dataset_id, gene_symbol, None, "sess", None)
    assert result["success"] == -1
    assert message in result["message"]
    assert validation.get_analysis_calls == []


def test_validate_args_returns_adata_and_analysis(validation):
    result = tsne_data.validate_args("DS1", "G1", {"id": "A1"}, "sess", None)
    assert result["success"] == 1
    assert result["adata"] is validation.adata
    assert result["ana"].dataset_id == "DS1"
    assert validation.get_analysis_calls == [{"analysis": {"id": "A1"}, "dataset_id": "DS1", "is_spatial": False}]
    assert validation.get_adata_calls == [{"backed": True}]


def test_validate_args_loads_spatial_datasets_in_memory(validation):
    assert tsne_data.validate_args("SP1", "G1", None, "sess", None)["success"] == 1
    assert validation.get_analysis_calls[0]["is_spatial"] is True
    assert validation.get_adata_calls == [{}]


@pytest.mark.parametrize("attr, message", [("analysis_error", "Analysis for this dataset is unavailable"),
                                           ("adata_error", "Could not create dataset object")])
def test_validate_args_reports_loader_failures(validation, attr, message):
    setattr(validation, attr, FileNotFoundError("gone"))
    result = tsne_data.validate_args("DS1", "G1", None, "sess", None)
    assert result == {"success": -1, "message": result["message"]}
    assert message in result["message"]


def test_validate_args_reports_projection_errors(validation, monkeypatch):
    def fail(adata, dataset_id, projection_id):
        raise PlotError("Projection P1 is not available")

    monkeypatch.setattr(tsne_data, "create_projection_adata", fail)
    result = tsne_data.validate_args("DS1", "G1", None, "sess", "P1")
    assert result == {"success": -1, "message": "Projection P1 is not available"}


def test_validate_args_requires_gene_symbol_column(validation):
    validation.adata.var = validation.adata.var.drop(columns="gene_symbol")
    result = tsne_data.validate_args("DS1", "G1", None, "sess", None)
    assert result["success"] == -1
    assert "gene_symbol column" in result["message"]


# ---------------------------------------------------------------- filter_adata_genes


@pytest.mark.parametrize(
    "genes, message",
    [(["NOPE"], 'The searched gene symbol "NOPE" could not be found'),
     (["NOPE", "NADA"], "The searched gene symbols could not be found")],
)
def test_filter_adata_genes_raises_when_no_gene_is_found(genes, message):
    with pytest.raises(ValueError, match=message):
        tsne_data.filter_adata_genes(build_adata(), genes)


def test_filter_adata_genes_keeps_requested_genes_as_dense_matrix():
    X = scipy.sparse.csr_matrix(np.arange(18, dtype=np.float32).reshape(6, 3))
    adata = build_adata(gene_symbols=("G1", "G2", "G3"), X=X)
    selected = tsne_data.filter_adata_genes(adata, ["G1", "G3", "NOPE"])
    assert selected.var.gene_symbol.tolist() == ["G1", "G3"]
    assert isinstance(selected.X, np.ndarray)
    assert np.array_equal(selected.X, X.toarray()[:, [0, 2]])
    assert list(selected.obs.columns) == ["cluster", "tSNE_1", "tSNE_2"]


def test_filter_adata_genes_closes_backed_file(tmp_path):
    path = tmp_path / "DS1.h5ad"
    build_adata().write_h5ad(path)
    adata = anndata.read_h5ad(path, backed="r")
    selected = tsne_data.filter_adata_genes(adata, ["G2"])
    assert not adata.file.is_open
    assert not selected.isbacked
    assert selected.shape == (6, 1)


# ---------------------------------------------------------------- validate_analysis_inputs


def test_obs_columns_become_embeddings_for_primary_dataset():
    adata = tsne_data.validate_analysis_inputs(build_adata(), "DS1", "DS1", "tSNE_1", "tSNE_2")
    expected = adata.obs[["tSNE_1", "tSNE_2"]].to_numpy()
    for key in ("X_tsne", "X_umap", "X_pca"):
        assert np.array_equal(adata.obsm[key], expected)


@pytest.mark.parametrize("analysis_id", [None, "null", "undefined"])
def test_unset_analysis_counts_as_primary(analysis_id):
    adata = tsne_data.validate_analysis_inputs(build_adata(), analysis_id, "DS1", "tSNE_1", "tSNE_2")
    assert "X_tsne" in adata.obsm


def test_missing_obs_column_is_rejected():
    with pytest.raises(ValueError, match="Dataseries UMAP_2 was selected but not present"):
        tsne_data.validate_analysis_inputs(build_adata(), "DS1", "DS1", "tSNE_1", "UMAP_2")


@pytest.mark.parametrize(
    "x_axis, y_axis, obsm_key, message",
    [("X_tsne_1", "X_tsne_2", "X_tsne", "tSNE"), ("X_umap_1", "X_umap_2", "X_umap", "UMAP"),
     ("X_pca_1", "X_pca_2", "X_pca", "PCA")],
)
def test_analysis_axes_come_from_obsm(x_axis, y_axis, obsm_key, message):
    with pytest.raises(ValueError, match=f"Analysis {message}.*not present in adata.obsm"):
        tsne_data.validate_analysis_inputs(build_adata(), "DS1", "DS1", x_axis, y_axis)

    adata = build_adata()
    embedding = np.ones((adata.n_obs, 2))
    adata.obsm[obsm_key] = embedding
    adata = tsne_data.validate_analysis_inputs(adata, "DS1", "DS1", x_axis, y_axis)
    assert adata.obsm[obsm_key] is embedding


def test_secondary_analysis_is_not_checked():
    adata = tsne_data.validate_analysis_inputs(build_adata(), "A1", "DS1", "X_tsne_1", "X_tsne_2")
    assert len(adata.obsm) == 0


# ---------------------------------------------------------------- prepare_adata


@pytest.mark.parametrize("flip_x, flip_y", [(True, False), (False, True), (True, True)])
def test_prepare_adata_flips_every_embedding(flip_x, flip_y):
    adata = build_adata()
    original = np.arange(adata.n_obs * 2, dtype=float).reshape(-1, 2)
    for key in ("X_tsne", "X_umap", "X_pca"):
        adata.obsm[key] = original.copy()
    adata = tsne_data.prepare_adata(adata, flip_x=flip_x, flip_y=flip_y)
    sign = np.array([-1 if flip_x else 1, -1 if flip_y else 1])
    for key in ("X_tsne", "X_umap", "X_pca"):
        assert np.array_equal(adata.obsm[key], original * sign)


def test_prepare_adata_reorders_categories_and_skips_numeric_columns():
    adata = tsne_data.prepare_adata(build_adata(), order={"cluster": ["c", "a", "b"], "tSNE_1": [1, 2]})
    assert adata.obs["cluster"].cat.categories.tolist() == ["c", "a", "b"]
    assert adata.obs["cluster"].cat.ordered


def test_prepare_adata_applies_every_filter():
    adata = build_adata()
    adata.obs["batch"] = pd.Categorical(["x", "y"] * 3)
    adata = tsne_data.prepare_adata(adata, filters={"cluster": ["a", "b"], "batch": ["x"]})
    assert adata.obs_names.tolist() == ["cell0", "cell4"]


@pytest.mark.parametrize("kwargs, message", [
    ({"order": {"cluster": ["c", "a"]}}, "Sort order for series cluster does not match the categories"),
    ({"order": {"cluster": ["c", "a", "b", "d"]}}, "Sort order for series cluster does not match the categories"),
    ({"order": {"celltype": ["a"]}}, "Sort order series celltype not found in observation metadata"),
    ({"filters": {"celltype": ["a"]}}, "Filter series celltype not found in observation metadata"),
])
def test_prepare_adata_reports_curation_that_no_longer_fits_the_dataset(kwargs, message):
    # These were a ValueError from reorder_categories or a KeyError from obs, so the endpoint gave a 500
    with pytest.raises(PlotError, match=message) as excinfo:
        tsne_data.prepare_adata(build_adata(), **kwargs)
    assert str(excinfo.value).endswith("Please update curation.")


# ---------------------------------------------------------------- modify_genes_found_as_obs_columns


def test_obs_column_named_like_a_gene_is_renamed():
    adata = build_adata()
    adata.obs["G1"] = list("uvwxyz")
    adata = tsne_data.modify_genes_found_as_obs_columns(adata, ["G1", "G2"])
    assert "G1" not in adata.obs
    assert adata.obs["G1_orig"].tolist() == list("uvwxyz")
    assert list(adata.obs.columns) == ["cluster", "tSNE_1", "tSNE_2", "G1_orig"]


# ---------------------------------------------------------------- dedup_genes


def test_dedup_genes_without_duplicates(tmp_path):
    result = tsne_data.dedup_genes(build_adata(), analysis(tmp_path))
    assert (result["success"], result["message"]) == (1, "")
    assert result["adata"].var_names.tolist() == ["G1", "G2"]
    assert result["adata"].var["ensembl_id"].tolist() == ["ENSG0", "ENSG1"]
    assert result["dedup_copy"] == tmp_path / "DS1.dups_removed.h5ad"
    assert not result["dedup_copy"].exists()


def test_dedup_genes_keeps_first_ensembl_id(tmp_path):
    adata = build_adata(gene_symbols=("G1", "G2", "G1"))
    result = tsne_data.dedup_genes(adata, analysis(tmp_path))
    assert result["success"] == 2
    assert "duplicated across multiple Ensembl IDs" in result["message"]
    assert result["adata"].var_names.tolist() == ["G1", "G2"]
    assert result["adata"].var["ensembl_id"].tolist() == ["ENSG0", "ENSG1"]


@pytest.mark.parametrize("dataset_name", ["DS1.zarr", "DS1.h5ad"])
def test_spatial_dedup_copy_never_reuses_the_dataset_path(tmp_path, dataset_name):
    dataset_path = tmp_path / dataset_name
    result = tsne_data.dedup_genes(build_adata(), NonPrimarySpatialAnalysis(dataset_path))
    assert result["dedup_copy"] == tmp_path / "DS1.dups_removed.h5ad"


def test_plotting_a_saved_spatial_analysis_keeps_its_h5ad(tmp_path):
    """The dedup cleanup used to delete the analysis' own .h5ad (non-primary spatial analyses aren't .zarr)."""
    dataset_path = tmp_path / "DS1.h5ad"
    dataset_path.write_bytes(b"analysis")
    adata = build_adata()
    adata.obsm["X_tsne"] = adata.obs[["tSNE_1", "tSNE_2"]].to_numpy()
    result = tsne_data.generate_tsne_figure(
        adata, NonPrimarySpatialAnalysis(dataset_path), ["G1"], "tsne_static", "X_tsne_1", "X_tsne_2"
    )
    assert result["success"] == 1, result
    assert dataset_path.read_bytes() == b"analysis"


# ---------------------------------------------------------------- TSNEData / MGTSNEData endpoints


@pytest.fixture
def client(validation):
    app = Flask(__name__)
    api = Api(app)
    api.add_resource(tsne_data.TSNEData, "/plot/<dataset_id>/tsne")
    api.add_resource(tsne_data.MGTSNEData, "/plot/<dataset_id>/mg_tsne")
    return app.test_client()


def post(client, route, **body):
    resp = client.post(f"/plot/DS1/{route}", json={"plot_type": "tsne_static", **body})
    assert resp.status_code == 200, resp.get_data(as_text=True)
    return resp.get_json()


@pytest.mark.parametrize(
    "route, body, message",
    [
        ("tsne", {}, "Request needs a gene symbol"),
        ("mg_tsne", {}, "Request needs a gene symbol"),
        ("tsne", {"gene_symbol": "NOPE"}, 'The searched gene symbol "NOPE" could not be found'),
        ("mg_tsne", {"gene_symbols": ["NOPE", "NADA"]}, "The searched gene symbols could not be found"),
        ("tsne", {"gene_symbol": "G1", "x_axis": "UMAP_1"}, "Dataseries UMAP_1 was selected"),
        ("tsne", {"gene_symbol": "G1", "plot_type": "violin"}, "violin was not a valid plot type"),
        ("tsne", {"gene_symbol": "G1", "vmax": 0.2, "vmin": 0.5}, "Maximum cap value must be greater"),
        ("tsne", {"gene_symbol": "G1", "vmax": "0.2", "vmin": "0.5"}, "Maximum cap value must be greater"),
        ("tsne", {"gene_symbol": "G1", "vmax": 0, "vmin": 0}, "Maximum cap value must be greater"),
        ("tsne", {"gene_symbol": "G1", "vmin": "abc"}, "The vmin and vmax parameters must be numbers."),
        ("mg_tsne", {"gene_symbols": ["G1", "G2"], "vmax": "x"}, "The vmin and vmax parameters must be numbers."),
        ("tsne", {"gene_symbol": "G1", "center_around_median": True, "vmax": 0.0001}, "Median expression value"),
        ("tsne", {"gene_symbol": "G1", "center_around_median": True, "vmin": 0.9999}, "Median expression value"),
        # Saved curations that no longer fit the dataset (these were 500s)
        ("tsne", {"gene_symbol": "G1", "order": {"cluster": ["c", "a"]}},
         "Sort order for series cluster does not match the categories in the dataset. Please update curation."),
        ("mg_tsne", {"gene_symbols": ["G1", "G2"], "order": {"celltype": ["a"]}},
         "Sort order series celltype not found in observation metadata for dataset. Please update curation."),
        ("tsne", {"gene_symbol": "G1", "obs_filters": {"celltype": ["a"]}},
         "Filter series celltype not found in observation metadata for dataset. Please update curation."),
        ("mg_tsne", {"gene_symbols": ["G1", "G2"], "obs_filters": {"celltype": ["a"]}},
         "Filter series celltype not found in observation metadata for dataset. Please update curation."),
    ],
)
def test_endpoint_errors(client, route, body, message):
    data = post(client, route, **body)
    assert data["success"] == -1
    assert message in data["message"]
    assert "image" not in data


def test_unknown_dataset(client):
    resp = client.post("/plot/missing/tsne", json={"gene_symbol": "G1"})
    assert resp.get_json() == {"success": -1, "message": "No dataset found with that ID"}


@pytest.mark.parametrize("route, body", [("tsne", {"gene_symbol": "G1"}),
                                         ("mg_tsne", {"gene_symbols": ["G1", "G2"], "colorize_legend_by": "cluster"})])
def test_endpoint_returns_image(client, route, body):
    data = post(client, route, **body)
    assert data["success"] == 1, data
    assert data["image_format"] == "webp"
    assert data["image"]
