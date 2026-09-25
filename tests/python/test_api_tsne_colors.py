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
