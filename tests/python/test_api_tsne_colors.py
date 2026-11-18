"""tsne_data.generate_tsne_figure: category colors passed to scanpy for static embeddings."""

import types

import anndata
import numpy as np
import pandas as pd
import pytest

from resources import tsne_data


class _StopBeforePlotting(Exception):
    pass


def build_adata(categories):
    n_cells = 3 * len(categories)
    rng = np.random.default_rng(0)
    obs = pd.DataFrame(
        {
            "cluster": pd.Categorical([c for c in categories for _ in range(3)], categories=categories),
            "tSNE_1": rng.normal(size=n_cells),
            "tSNE_2": rng.normal(size=n_cells),
        },
        index=[f"cell{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame({"gene_symbol": ["G1", "G2"]}, index=["ENSG1", "ENSG2"])
    return anndata.AnnData(X=rng.uniform(size=(n_cells, 2)).astype(np.float32), obs=obs, var=var)


@pytest.fixture
def plotted_uns(monkeypatch):
    """Run generate_tsne_figure and return the .uns scanpy would have plotted with."""
    captured = {}

    def fake_embedding(adata, **kwargs):
        captured.update(adata.uns)
        raise _StopBeforePlotting

    monkeypatch.setattr(tsne_data.sc.pl, "embedding", fake_embedding)

    def run(adata, **kwargs):
        ana = types.SimpleNamespace(dataset_id="DS1", id="DS1", dataset_path="/nonexistent/DS1.h5ad")
        with pytest.raises(_StopBeforePlotting):
            tsne_data.generate_tsne_figure(
                adata, ana, ["G1"], "tsne_static", "tSNE_1", "tSNE_2", colorize_by="cluster", **kwargs
            )
        return captured

    return run


@pytest.mark.parametrize("categories", [["a"], ["a", "b"], ["a", "b", "c"]])
def test_curator_colors_used_for_any_number_of_categories(plotted_uns, categories):
    colors = {c: f"#0000{i:02x}" for i, c in enumerate(categories)}
    uns = plotted_uns(build_adata(categories), colors=colors)
    assert uns["cluster_colors"] == [colors[c] for c in categories]


def test_incomplete_curator_colors_fall_back_to_obs_color_column(plotted_uns):
    adata = build_adata(["a", "b"])
    adata.obs["cluster_colors"] = ["#111111"] * 3 + ["#222222"] * 3
    uns = plotted_uns(adata, colors={"a": "#ff0000"})
    assert uns["cluster_colors"] == ["#111111", "#222222"]


def test_colorblind_mode_with_a_single_category(plotted_uns):
    uns = plotted_uns(build_adata(["a"]), colorblind_mode=True)
    assert len(uns["cluster_colors"]) == 1


@pytest.mark.parametrize("colorblind_mode, palette", [(False, "YlOrRd"), (True, "cividis_r")])
def test_zero_gray_expression_scale_keeps_the_chosen_palette(monkeypatch, colorblind_mode, palette):
    """make_zero_gray (on by default) must not replace the colorblind cividis scale."""
    import matplotlib.pyplot as plt

    captured = {}

    def fake_embedding(adata, **kwargs):
        captured.update(kwargs)
        raise _StopBeforePlotting

    monkeypatch.setattr(tsne_data.sc.pl, "embedding", fake_embedding)
    ana = types.SimpleNamespace(dataset_id="DS1", id="DS1", dataset_path="/nonexistent/DS1.h5ad")
    with pytest.raises(_StopBeforePlotting):
        tsne_data.generate_tsne_figure(
            build_adata(["a"]), ana, ["G1"], "tsne_static", "tSNE_1", "tSNE_2",
            expression_palette="YlOrRd", colorblind_mode=colorblind_mode, make_zero_gray=True,
        )

    cmap = captured["color_map"]
    assert np.allclose(cmap(0.0), [192 / 256, 192 / 256, 192 / 256, 1])
    assert np.allclose(cmap(1.0), plt.get_cmap(palette)(1.0))


def build_adata_with_score(categories=("a", "b")):
    """build_adata plus a continuous obs column whose range (100-200) is far from expression (0-1)."""
    adata = build_adata(list(categories))
    adata.obs["score"] = np.linspace(100, 200, adata.n_obs)
    return adata


def captured_embedding_kwargs(monkeypatch, adata, **kwargs):
    """Run generate_tsne_figure and return the kwargs it passed to scanpy."""
    captured = {}

    def fake_embedding(adata, **kwargs):
        captured.update(kwargs)
        raise _StopBeforePlotting

    monkeypatch.setattr(tsne_data.sc.pl, "embedding", fake_embedding)
    ana = types.SimpleNamespace(dataset_id="DS1", id="DS1", dataset_path="/nonexistent/DS1.h5ad")
    with pytest.raises(_StopBeforePlotting):
        tsne_data.generate_tsne_figure(adata, ana, ["G1"], "tsne_static", "tSNE_1", "tSNE_2", **kwargs)
    return captured


def test_expression_caps_skip_continuous_colorize_panel(monkeypatch):
    kwargs = captured_embedding_kwargs(
        monkeypatch, build_adata_with_score(), colorize_by="score", vmin=0.1, vmax=0.5
    )
    assert kwargs["color"] == ["G1", "score"]
    assert kwargs["vmin"] == [0.1, None]
    assert kwargs["vmax"] == [0.5, None]
    assert kwargs["vcenter"] == [None, None]


def test_expression_caps_unchanged_for_categorical_colorize(monkeypatch):
    kwargs = captured_embedding_kwargs(
        monkeypatch, build_adata_with_score(), colorize_by="cluster", vmin=0.1, vmax=0.5
    )
    assert kwargs["vmin"] == 0.1
    assert kwargs["vmax"] == 0.5


@pytest.fixture
def plotted_figure(monkeypatch):
    """Run generate_tsne_figure with the real scanpy plot and return the figure it drew."""
    real_embedding = tsne_data.sc.pl.embedding
    captured = {}

    def keep_figure(adata, **kwargs):
        captured["fig"] = real_embedding(adata, **kwargs)
        return captured["fig"]

    monkeypatch.setattr(tsne_data.sc.pl, "embedding", keep_figure)

    def run(adata, **kwargs):
        ana = types.SimpleNamespace(dataset_id="DS1", id="DS1", dataset_path="/nonexistent/DS1.h5ad")
        result = tsne_data.generate_tsne_figure(
            adata, ana, ["G1"], "tsne_static", "tSNE_1", "tSNE_2", **kwargs
        )
        assert result["success"] == 1, result
        return captured["fig"]

    return run


def panel_scatter(fig, title):
    """Return the scatter of the panel with the given title."""
    from matplotlib.collections import PathCollection

    axes = next(a for a in fig.get_axes() if a.get_title() == title)
    return next(c for c in axes.collections if isinstance(c, PathCollection))


@pytest.mark.parametrize(
    "expression_palette, colorblind_mode, expected",
    [("YlOrRd", False, "viridis"), ("viridis", False, "plasma"), ("viridis", True, "viridis")],
)
def test_continuous_colorize_panel_has_its_own_colormap(
    plotted_figure, expression_palette, colorblind_mode, expected
):
    fig = plotted_figure(
        build_adata_with_score(),
        colorize_by="score",
        expression_palette=expression_palette,
        colorblind_mode=colorblind_mode,
        vmax=0.5,
    )
    gene = panel_scatter(fig, "G1")
    score = panel_scatter(fig, "score")

    assert score.get_cmap().name == expected
    assert score.colorbar.cmap.name == expected
    assert gene.get_cmap().name != expected
    # Expression keeps its cap; the annotation uses its own range
    assert gene.norm.vmax == pytest.approx(0.5)
    assert score.norm.vmax == pytest.approx(200)


def test_plot_by_group_with_continuous_colorize(plotted_figure):
    """Group subplots use the expression colormap and range; the colorize_by panel keeps its own."""
    adata = build_adata_with_score()
    max_expression = float(adata[:, "ENSG1"].X.max())
    fig = plotted_figure(adata, colorize_by="score", plot_by_group="cluster", expression_palette="YlOrRd")

    titles = [a.get_title() for a in fig.get_axes() if a.get_label() != "<colorbar>"]
    assert titles == ["a", "b", "G1", "score"]

    for title in ["a", "b", "G1"]:
        assert panel_scatter(fig, title).get_cmap().name != "viridis"
    for title in ["a", "b"]:
        assert panel_scatter(fig, title).norm.vmax == pytest.approx(max_expression)

    score = panel_scatter(fig, "score")
    assert score.get_cmap().name == "viridis"
    assert score.norm.vmax == pytest.approx(200)
