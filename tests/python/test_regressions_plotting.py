"""Regression tests for closed GitHub issues in gear.plotting.generate_plot (single-gene Plotly displays)."""

import pandas as pd
import pytest

pytest.importorskip("plotly")

from gear import plotting  # noqa: E402


def plot(df, plot_type, **kwargs):
    kwargs = {"x_title": "X", "y_title": "Y", **kwargs}
    return plotting.generate_plot(df.copy(), "x", "y", plot_type=plot_type, **kwargs)


def x_axis_types(fig) -> set:
    return {getattr(fig.layout, key).type for key in fig.layout if key.startswith("xaxis")}


# ---------------------------------------------------------------- #1293 bar charts with replicate rows


def test_bar_plot_with_replicate_rows_has_one_mean_bar_per_group():
    """Regression for #1293: replicate observations (no "replicate" column) give one bar per x/color group at the
    mean, side by side, instead of one stacked segment per observation."""
    df = pd.DataFrame({
        "x": pd.Categorical(["ctrl"] * 4 + ["treated"] * 4),
        "sex": pd.Categorical(["Female", "Female", "Male", "Male"] * 2),
        "y": [1.0, 3.0, 2.0, 4.0, 5.0, 7.0, 6.0, 10.0],
    })
    fig = plot(df, "bar", color_name="sex", colormap={"Female": "#ec8ec9", "Male": "#0fb7ff"})
    bars = {t.name: dict(zip(t.x, t.y)) for t in fig.data}
    assert bars == {"Female": {"ctrl": 2.0, "treated": 6.0}, "Male": {"ctrl": 3.0, "treated": 8.0}}
    assert all(len(t.x) == len(set(t.x)) for t in fig.data)
    assert fig.layout.barmode == "group"


# ---------------------------------------------------------------- #1285 line plots with numeric x


@pytest.mark.parametrize("dtype", ["float64", "int64", "float32", "int32", "Int64", "Float64"])
def test_line_plot_numeric_x_is_sorted(dtype):
    """Regression for #1285: a line plot with a numeric x axis connects points in ascending x order (0.0 last made the line loop back)."""
    df = pd.DataFrame({"x": pd.Series([2, 0, 3, 1], dtype=dtype), "y": [20.0, 0.0, 30.0, 10.0]})
    fig = plot(df, "line")
    assert [list(t.x) for t in fig.data] == [[0, 1, 2, 3]]
    assert list(fig.data[0].y) == [0.0, 10.0, 20.0, 30.0]


def test_line_plot_numeric_x_is_sorted_per_color_group():
    """Regression for #1285: with a color group, each line's x values are ascending too."""
    df = pd.DataFrame({
        "x": [1.0, 0.0, 2.0, 2.0, 0.0, 1.0],
        "g": pd.Categorical(["p", "p", "p", "q", "q", "q"]),
        "y": [1.0, 0.0, 2.0, 5.0, 3.0, 4.0],
    })
    fig = plot(df, "line", color_name="g")
    assert {t.name: (list(t.x), list(t.y)) for t in fig.data} == {
        "p": ([0.0, 1.0, 2.0], [0.0, 1.0, 2.0]),
        "q": ([0.0, 1.0, 2.0], [3.0, 4.0, 5.0]),
    }


# ---------------------------------------------------------------- #1071 categorical x axis


@pytest.mark.parametrize("plot_type", ["violin", "bar", "box", "strip"])
@pytest.mark.parametrize("categories", [[0, 1, 2], ["0", "a", "b"]], ids=["ints", "zero-first"])
def test_categorical_plot_types_force_category_x_axis(plot_type, categories):
    """Regression for #1071: violin, bar, box and strip plots use a categorical x axis even when the first x value is a number."""
    df = pd.DataFrame({
        "x": pd.Categorical([c for c in categories for _ in range(2)]),
        "g": pd.Categorical(["p", "q"] * len(categories)),
        "y": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    })
    fig = plot(df, plot_type, color_name="g")
    assert x_axis_types(fig) == {"category"}
    assert {str(v) for t in fig.data for v in t.x} == {str(c) for c in categories}


def test_violin_with_integer_x_and_no_other_grouping():
    """Regression for #1071: a violin plot of integer x categories with no color or facet grouping renders on a categorical axis."""
    df = pd.DataFrame({"x": pd.Categorical([0, 0, 1, 1]), "y": [1.0, 2.0, 3.0, 4.0]})
    fig = plot(df, "violin")
    assert x_axis_types(fig) == {"category"}


# ---------------------------------------------------------------- #1145 violin facet rows with an empty group


@pytest.mark.parametrize("color_name", [None, "g"])
def test_violin_facet_rows_get_placeholder_for_empty_group(color_name):
    """Regression for #1145: when a facet row has no cells for an x group, a hidden placeholder trace keeps
    that row's violins aligned with the other rows. (test_plotting_helpers covers the unfaceted case.)"""
    df = pd.DataFrame({
        "r": pd.Categorical(["r1"] * 4 + ["r2"] * 2),
        "x": pd.Categorical(["a", "a", "b", "b", "a", "a"]),
        "g": pd.Categorical(["p", "q"] * 3),
        "y": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    })
    fig = plot(df, "violin", facet_row="r", color_name=color_name)

    rows = {"y": "r1", "y2": "r2"}
    per_row = {}
    for t in fig.data:
        per_row.setdefault(rows[t.yaxis], []).append((t.name, tuple(t.x), t.opacity))

    colors = ["p", "q"] if color_name else [None]
    # every facet row has a trace for every (color, x) combination, in the same order
    assert [(n, x[0]) for n, x, _ in per_row["r1"]] == [(c, x) for c in colors for x in "ab"]
    assert [(n, x[0]) for n, x, _ in per_row["r2"]] == [(c, x) for c in colors for x in "ab"]
    # only r2's "b" traces are placeholders, and they are invisible
    hidden = [(n, x) for n, x, opacity in per_row["r1"] + per_row["r2"] if opacity == 0.0]
    assert hidden == [(c, ("b",)) for c in colors]
