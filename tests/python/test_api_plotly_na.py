"""POST /plot/<dataset_id>: missing obs values are plotted as an "NA" group (issue #888)."""

import types

import anndata
import numpy as np
import pandas as pd
import pytest
from flask import Flask
from flask_restful import Api

from resources import plotly_data

HAIR, SUPPORT = "#1f77b4", "#ff7f0e"


def build_obs():
    """9 cells: 3 hair, 3 support and 3 with no cluster (so no cluster color either)."""
    clusters = ["hair"] * 3 + ["support"] * 3 + [None] * 3
    colors = {"hair": HAIR, "support": SUPPORT}
    return pd.DataFrame(
        {
            "cluster": pd.Categorical(clusters, categories=["support", "hair"], ordered=True),
            "cluster_colors": pd.Categorical([colors.get(c) for c in clusters]),
            "batch": pd.Categorical(["b1", "b2", "b3"] * 3),      # no missing values
            "n_genes": [100.0, np.nan] + [200.0] * 7,             # numeric: left as NaN
            "tSNE_1": np.arange(9, dtype=float),
            "tSNE_2": np.arange(9, dtype=float),
        },
        index=[f"cell{i}" for i in range(9)],
    )


def build_adata():
    var = pd.DataFrame({"gene_symbol": ["Myo7a", "Gapdh"]}, index=["ENSMUSG00000030761", "ENSMUSG00000057666"])
    X = np.arange(18, dtype=np.float32).reshape(9, 2)
    return anndata.AnnData(X=X, obs=build_obs(), var=var)


class TestFillMissingObs:
    def test_categorical_gets_na_category_last_and_keeps_order(self):
        obs = plotly_data.fill_missing_obs(build_obs())
        assert obs["cluster"].tolist() == ["hair"] * 3 + ["support"] * 3 + ["NA"] * 3
        assert obs["cluster"].cat.categories.tolist() == ["support", "hair", "NA"]
        assert obs["cluster"].cat.ordered

    def test_color_column_gets_missing_value_color(self):
        obs = plotly_data.fill_missing_obs(build_obs())
        assert obs["cluster_colors"].tolist() == [HAIR] * 3 + [SUPPORT] * 3 + [plotly_data.NA_COLOR] * 3

    def test_complete_and_numeric_columns_unchanged(self):
        obs = plotly_data.fill_missing_obs(build_obs())
        assert "NA" not in obs["batch"].cat.categories
        assert np.isnan(obs["n_genes"].iloc[1])

    def test_object_column_filled(self):
        obs = pd.DataFrame({"label": pd.Series(["x", None, "y"], dtype=object)})
        assert plotly_data.fill_missing_obs(obs)["label"].tolist() == ["x", "NA", "y"]

    def test_existing_na_category_reused(self):
        obs = pd.DataFrame({"c": pd.Categorical(["x", None], categories=["x", "NA"])})
        filled = plotly_data.fill_missing_obs(obs)["c"]
        assert filled.tolist() == ["x", "NA"]
        assert filled.cat.categories.tolist() == ["x", "NA"]


@pytest.fixture
def client(monkeypatch):
    mod = plotly_data
    monkeypatch.setattr(
        mod.geardb, "get_dataset_by_id",
        lambda dataset_id, *a, **k: types.SimpleNamespace(id=dataset_id, dtype="single-cell-rnaseq"),
    )
    monkeypatch.setattr(mod, "get_adata_from_analysis", lambda *a, **k: build_adata())

    app = Flask(__name__)
    Api(app).add_resource(mod.PlotlyData, "/plot/<dataset_id>")
    return app.test_client()


def plot(client, **body):
    body = {"gene_symbol": "Myo7a", "plot_type": "scatter", "y_axis": "raw_value", **body}
    resp = client.post("/plot/DS1", json=body)
    assert resp.status_code == 200, resp.get_data(as_text=True)
    result = resp.get_json()
    assert result["success"] == 1, result["message"]
    return result


def traces(result):
    return result["plot_json"]["data"]


def x_values(result):
    return [x for t in traces(result) for x in t["x"]]


def trace_colors(result):
    """legend group name -> marker color, for a discrete color scatter."""
    return {t["name"]: t["marker"]["color"] for t in traces(result)}


def test_na_group_plotted_on_x_axis(client):
    xs = x_values(plot(client, x_axis="cluster"))
    assert len(xs) == 9
    assert xs.count("NA") == 3


def test_na_group_gets_missing_value_color_from_color_column(client):
    colors = trace_colors(plot(client, x_axis="tSNE_1", y_axis="tSNE_2", color_name="cluster"))
    assert colors == {"hair": HAIR, "support": SUPPORT, "NA": plotly_data.NA_COLOR}


def test_curator_colors_without_na_are_not_shifted(client):
    """A saved color map without an "NA" entry keeps its colors, and NA gets the missing-value color."""
    colors = trace_colors(plot(
        client, x_axis="tSNE_1", y_axis="tSNE_2", color_name="cluster",
        colors={"hair": "#00ff00", "support": "#0000ff"},
    ))
    assert colors == {"hair": "#00ff00", "support": "#0000ff", "NA": plotly_data.NA_COLOR}


def test_curator_color_for_na_is_used(client):
    colors = trace_colors(plot(
        client, x_axis="tSNE_1", y_axis="tSNE_2", color_name="cluster",
        colors={"hair": "#00ff00", "support": "#0000ff", "NA": "#ff00ff"},
    ))
    assert colors["NA"] == "#ff00ff"


def test_filter_on_na(client):
    xs = x_values(plot(client, x_axis="cluster", obs_filters={"cluster": ["hair", "NA"]}))
    assert sorted(set(xs)) == ["NA", "hair"]
    assert len(xs) == 6


@pytest.mark.parametrize("order", [["NA", "hair", "support"], ["hair", "support"]])
def test_saved_order_with_or_without_na(client, order):
    """The curator UI lists "NA" as a level, so saved orders may or may not include it."""
    result = plot(client, x_axis="cluster", order={"cluster": order})
    expected = order if "NA" in order else order + ["NA"]
    assert result["plot_order"]["cluster"] == expected
