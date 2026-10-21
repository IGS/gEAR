"""Helpers in gear.plotting (single-gene Plotly displays) and generate_plot smoke tests."""

import pandas as pd
import pytest

pytest.importorskip("plotly")

from gear import plotting  # noqa: E402

PURPLE = "#401362"


# ---------------------------------------------------------------- grouping and aggregation

@pytest.mark.parametrize("kwargs, expected", [
    ({}, []),
    ({"x": "x"}, ["x"]),
    ({"facet_row": "r", "facet_col": "c", "color_name": "g", "x": "x"}, ["r", "c", "g", "x"]),
    ({"facet_row": "r", "facet_col": "r", "x": "x"}, ["r", "x"]),
    ({"facet_col": "c", "x": "x"}, ["c", "x"]),
    ({"facet_row": "r", "color_name": "r", "x": "r"}, ["r"]),
    ({"x": "x", "y": "y"}, ["x", "y"]),
])
def test_build_priority_groups(kwargs, expected):
    assert plotting._build_priority_groups(**kwargs) == expected


@pytest.fixture
def expr_df():
    return pd.DataFrame({
        "x": pd.Categorical(["a", "a", "b", "b"]),
        "g": pd.Categorical(["p", "q", "p", "p"]),
        "y": [1.0, 3.0, 5.0, 5.0],
    })


def test_aggregate_dataframe_mean_and_std(expr_df):
    agg = plotting._aggregate_dataframe(expr_df, "x", "y")
    assert agg.columns.tolist() == ["x", "y", "std"]
    assert agg["y"].tolist() == [2.0, 5.0]
    assert agg["std"].tolist() == pytest.approx([1.4142, 0.0], abs=1e-4)


def test_aggregate_dataframe_by_color_drops_empty_groups_and_zero_fills_std(expr_df):
    agg = plotting._aggregate_dataframe(expr_df, "x", "y", color_name="g")
    assert list(zip(agg["g"], agg["x"], agg["y"], agg["std"])) == [("p", "a", 1.0, 0.0), ("p", "b", 5.0, 0.0), ("q", "a", 3.0, 0.0)]


def test_aggregate_dataframe_without_groups_is_unchanged(expr_df):
    assert plotting._aggregate_dataframe(expr_df, None, "y") is expr_df


def test_aggregate_dataframe_continuous_color_is_averaged_per_group():
    df = pd.DataFrame({"x": pd.Categorical(["a", "a", "b", "b"]), "y": [1.0, 3.0, 5.0, 5.0], "c": [1.0, 2.0, 3.0, 4.0]})
    agg = plotting._aggregate_dataframe(df, "x", "y", color_name="c")
    assert list(agg.columns) == ["x", "c", "y", "std"]
    assert agg["x"].tolist() == ["a", "b"]
    assert agg["c"].tolist() == [1.5, 3.5]
    assert agg["y"].tolist() == [2.0, 5.0]
    assert agg["std"].tolist() == pytest.approx([2 ** 0.5, 0.0])


def test_aggregate_dataframe_continuous_color_same_as_x():
    df = pd.DataFrame({"x": [1.0, 1.0, 2.0], "y": [1.0, 3.0, 5.0]})
    agg = plotting._aggregate_dataframe(df, "x", "y", color_name="x")
    assert agg["x"].tolist() == [1.0, 2.0]
    assert agg["y"].tolist() == [2.0, 5.0]


# ---------------------------------------------------------------- colors

def test_adjust_colorscale_dict_is_discrete_map():
    args = plotting._adjust_colorscale({}, {"a": "red"})
    assert args == {"color_discrete_map": {"a": "red"}}


def test_adjust_colorscale_list_is_continuous_scale():
    scale = [[0, "white"], [1, "black"]]
    assert plotting._adjust_colorscale({}, scale) == {"color_continuous_scale": scale}


def test_adjust_colorscale_defaults_to_purples():
    scale = plotting._adjust_colorscale({})["color_continuous_scale"]
    assert scale[0] == [0, "rgb(218, 183, 193)"] and scale[-1] == [1, "rgb(63, 19, 98)"]


@pytest.mark.parametrize("palette, first", [
    ("viridis", None),
    ("multicolor_diverging", [0.0, "rgb(154, 14, 234)"]),
    ("bublrd", [0, "rgb(173, 216, 230)"]),
])
def test_adjust_colorscale_palette(palette, first):
    scale = plotting._adjust_colorscale({}, palette=palette)["color_continuous_scale"]
    if first is None:
        assert scale == palette
    else:
        assert scale[0] == first


def test_plotly_color_map_single_group_is_gear_purple():
    assert plotting.plotly_color_map(["only"]) == {"only": PURPLE}


def test_plotly_color_map_cycles_d3_category20():
    names = [f"g{i}" for i in range(22)]
    colors = plotting.plotly_color_map(names)
    assert len(set(colors.values())) == 20
    assert colors["g0"] == colors["g20"] == "#1f77b4" and colors["g1"] == colors["g21"]


@pytest.mark.parametrize("rgb, expected", [((64, 19, 98), PURPLE), (("255", "0", "16"), "#ff0010"), ((0, 0, 0), "#000000")])
def test_rgb_to_hex(rgb, expected):
    assert plotting.rgb_to_hex(*rgb) == expected


# ---------------------------------------------------------------- layout helpers

def test_annotation_shift_uses_longest_category():
    assert plotting._determine_annotation_shift(pd.Series(["ab", "abcdef"], dtype="category")) == -(6 * 3 + 30)


def test_annotation_shift_non_categorical_is_constant():
    assert plotting._determine_annotation_shift(pd.Series([1.0, 200000.0])) == -30


def test_annotation_shift_non_string_categories():
    assert plotting._determine_annotation_shift(pd.Series([1, 100], dtype="category")) == -(3 * 3 + 30)


def test_translate_and_scale_maps_onto_0_to_9():
    series = pd.Series([2.0, 4.0, 6.0])
    assert [plotting._translate_and_scale(series, v) for v in series] == [0.0, 4.5, 9.0]


def test_translate_and_scale_constant_series():
    series = pd.Series([5, 5, 5])
    assert [plotting._translate_and_scale(series, v) for v in series] == [0, 0, 0]


def test_truncate_ticktext():
    ticks, mapping = plotting._truncate_ticktext(["short", "averyverylongname", "averyverylongother"])
    assert ticks == ["short", "averyve...", "averyve...~1"]
    assert mapping == {"averyve...": "averyverylongname", "averyve...~1": "averyverylongother"}
    assert plotting._truncate_ticktext([]) == (None, {})


def test_truncate_ticktext_non_string_values():
    assert plotting._truncate_ticktext([1, 12345678901]) == (["1", "1234567..."], {"1234567...": "12345678901"})


# ---------------------------------------------------------------- generate_plot

def plot(df, plot_type, **kwargs):
    kwargs = {"x_title": "X", "y_title": "Y", **kwargs}
    return plotting.generate_plot(df.copy(), "x", "y", plot_type=plot_type, **kwargs)


def test_bar_plot_with_continuous_color(expr_df):
    df = expr_df.assign(c=range(len(expr_df)))
    fig = plot(df, "bar", color_name="c")
    assert [t.type for t in fig.data] == ["bar"]
    assert list(fig.data[0].y) == [2.0, 5.0]
    assert fig.layout.coloraxis.colorscale


def test_scatter_with_continuous_color_is_not_aggregated(expr_df):
    """plotly_data drops "replicate" for scatter, so each point keeps its own continuous color."""
    df = expr_df.assign(c=[1.0, 2.0, 3.0, 4.0])
    fig = plot(df, "scatter", color_name="c")
    assert list(fig.data[0].y) == [1.0, 3.0, 5.0, 5.0]
    assert list(fig.data[0].marker.color) == [1.0, 2.0, 3.0, 4.0]


@pytest.mark.parametrize("bounds, expected", [
    ({"vmin": 1.5, "vmax": 3.5}, (1.5, 3.5)),
    ({"vmax": 3.0}, (1.0, 3.0)),     # the missing end comes from the data
    ({"vmin": 2.0}, (2.0, 4.0)),
])
def test_scatter_continuous_color_bounds(expr_df, bounds, expected):
    """#1166: vmin/vmax bound the colorbar of a continuous color series."""
    df = expr_df.assign(c=[1.0, 2.0, 3.0, 4.0])
    fig = plot(df, "scatter", color_name="c", palette="Viridis", **bounds)
    assert (fig.layout.coloraxis.cmin, fig.layout.coloraxis.cmax) == expected
    # The data itself is not clipped
    assert list(fig.data[0].marker.color) == [1.0, 2.0, 3.0, 4.0]


def test_scatter_without_bounds_leaves_coloraxis_auto(expr_df):
    df = expr_df.assign(c=[1.0, 2.0, 3.0, 4.0])
    fig = plot(df, "scatter", color_name="c", palette="Viridis")
    assert fig.layout.coloraxis.cmin is None and fig.layout.coloraxis.cmax is None


def test_color_bounds_ignored_for_categorical_color(expr_df):
    fig = plot(expr_df, "scatter", color_name="g", colormap={"p": "red", "q": "blue"}, vmin=0, vmax=1)
    assert fig.layout.coloraxis.cmin is None


def test_color_bound_beyond_data_is_plot_error(expr_df):
    df = expr_df.assign(c=[1.0, 2.0, 3.0, 4.0])
    with pytest.raises(plotting.PlotError, match="must be less than the maximum"):
        plot(df, "scatter", color_name="c", palette="Viridis", vmin=5.0)


def test_contour_color_bounds():
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0, 4.0], "y": [1.0, 2.0, 1.0, 2.0], "raw_value": [0.0, 1.0, 2.0, 8.0]})
    fig = plot(df, "contour", z="raw_value", vmax=4.0)
    assert (fig.data[0].zmin, fig.data[0].zmax) == (0.0, 4.0)


def test_bar_plot_aggregates_with_error_bars(expr_df):
    fig = plot(expr_df, "bar")
    assert [t.type for t in fig.data] == ["bar"]
    assert list(fig.data[0].y) == [2.0, 5.0]
    assert fig.data[0].error_y.array is not None
    assert fig.data[0].marker.color == PURPLE


def test_bar_plot_colored_by_group_has_trace_per_group(expr_df):
    fig = plot(expr_df, "bar", color_name="g", colormap={"p": "red", "q": "blue"})
    assert sorted((t.name, t.marker.color) for t in fig.data) == [("p", "red"), ("q", "blue")]


def test_scatter_plot(expr_df):
    fig = plot(expr_df, "scatter")
    assert [t.type for t in fig.data] == ["scatter"]
    assert len(fig.data[0].x) == 4


def test_scatter_plot_sized_by_constant_column():
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0], "s": [2.0, 2.0, 2.0]})
    fig = plot(df, "scatter", traces={"marker": {"size": "s", "sizemin": 3}})
    assert list(fig.data[0].marker.size) == [3, 3, 3]


def test_violin_plot_trace_per_group(expr_df):
    fig = plot(expr_df, "violin", color_name="g")
    assert {t.type for t in fig.data} == {"violin"}
    # one trace per (color, x) combination, the combination without data hidden
    assert len(fig.data) == 4
    assert sum(1 for t in fig.data if t.opacity == 0.0) == 1


def test_invalid_plot_type():
    with pytest.raises(plotting.PlotError, match="invalid"):
        plot(pd.DataFrame({"x": ["a"], "y": [1.0]}), "pie")


@pytest.mark.parametrize("kwargs", [{"facet_col": "f"}, {"category_orders": {"x": [1, 2]}}])
def test_bar_plot_with_numeric_categories(kwargs):
    df = pd.DataFrame({"x": pd.Categorical([1, 2, 1, 2]), "y": [1.0, 2.0, 3.0, 4.0], "f": pd.Categorical(["a", "a", "b", "b"])})
    fig = plot(df, "bar", **kwargs)
    assert {t.type for t in fig.data} == {"bar"}
