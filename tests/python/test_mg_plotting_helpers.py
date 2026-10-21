"""Pure helpers in gear.mg_plotting: condition validation, gene masks, axis labels, heatmap domains and palettes."""

import importlib
import sys
import types

import pandas as pd
import pytest


@pytest.fixture(scope="module")
def mg():
    """Import gear.mg_plotting, standing in for diffxpy (only used for differential expression) if it is missing."""
    pytest.importorskip("plotly")
    with pytest.MonkeyPatch.context() as mp:
        try:
            import diffxpy.api  # noqa: F401
        except ImportError:
            diffxpy = types.ModuleType("diffxpy")
            diffxpy.api = types.ModuleType("diffxpy.api")
            mp.setitem(sys.modules, "diffxpy", diffxpy)
            mp.setitem(sys.modules, "diffxpy.api", diffxpy.api)
        mp.delitem(sys.modules, "gear.mg_plotting", raising=False)
        yield importlib.import_module("gear.mg_plotting")


@pytest.fixture
def obs():
    return pd.DataFrame({"cell_type": ["a", "b", "c"], "batch": ["x", "y", "z"]})


# ---------------------------------------------------------------- quadrant / volcano conditions

def test_quadrant_conditions_split_into_key_and_values(mg, obs):
    result = mg.validate_quadrant_conditions(obs, "cell_type;-;a", "cell_type;-;b", "cell_type;-;c")
    assert result == ("cell_type", "a", "b", "c")


@pytest.mark.parametrize("conditions", [
    ("", "cell_type;-;b", "cell_type;-;c"),
    ("cell_type;-;a", None, "cell_type;-;c"),
    ("cell_type;-;a", "cell_type;-;b", ""),
])
def test_quadrant_conditions_require_all_three(mg, obs, conditions):
    with pytest.raises(mg.PlotError, match="three conditions .* quadrant plot"):
        mg.validate_quadrant_conditions(obs, *conditions)


@pytest.mark.parametrize("compare1, compare2", [
    ("batch;-;x", "cell_type;-;c"),
    ("cell_type;-;b", "batch;-;y"),
    ("batch;-;x", "batch;-;y"),
])
def test_quadrant_conditions_must_share_one_obs_column(mg, obs, compare1, compare2):
    with pytest.raises(mg.PlotError, match="same observation group"):
        mg.validate_quadrant_conditions(obs, "cell_type;-;a", compare1, compare2)


@pytest.mark.parametrize("bad", ["cell_type", "cell_type:a"])
def test_quadrant_condition_without_separator_is_a_plot_error(mg, obs, bad):
    with pytest.raises(mg.PlotError, match=f"'{bad}'"):
        mg.validate_quadrant_conditions(obs, "cell_type;-;a", bad, "cell_type;-;c")


def test_quadrant_condition_key_must_be_in_obs(mg, obs):
    with pytest.raises(mg.PlotError, match="missing not found"):
        mg.validate_quadrant_conditions(obs, "missing;-;a", "missing;-;b", "missing;-;c")


def test_volcano_conditions_split_into_key_and_values(mg, obs):
    assert mg.validate_volcano_conditions(obs, "cell_type;-;a", "cell_type;-;b") == ("cell_type", "a", "b")


def test_volcano_union_reference_is_shortened_to_rest(mg, obs):
    _, _, ref_val = mg.validate_volcano_conditions(obs, "cell_type;-;a", "cell_type;-;Union of the rest of the groups")
    assert ref_val == "rest"


@pytest.mark.parametrize("query, ref, match", [
    ("", "cell_type;-;b", "two conditions"),
    ("cell_type;-;a", "batch;-;x", "same observation group"),
    ("missing;-;a", "missing;-;b", "missing not found"),
    ("cell_type", "cell_type;-;b", "'cell_type'"),
    ("cell_type;-;a", "cell_type/b", "'cell_type/b'"),
])
def test_volcano_conditions_errors(mg, obs, query, ref, match):
    with pytest.raises(mg.PlotError, match=match):
        mg.validate_volcano_conditions(obs, query, ref)


# ---------------------------------------------------------------- gene searching

@pytest.fixture
def genes_df():
    return pd.DataFrame(
        {"gene_symbol": ["Gene1", "Gene2", "Dup", "Dup", "Other"]},
        index=["E1", "E2", "E3", "E4", "E5"],
    )


def test_gene_mask_requires_gene_symbol_column(mg):
    with pytest.raises(mg.PlotError, match="gene_symbol"):
        mg.create_dataframe_gene_mask(pd.DataFrame({"x": [1]}), ["a"])


def test_gene_mask_without_search_is_empty(mg, genes_df):
    assert mg.create_dataframe_gene_mask(genes_df, []) == (None, 1, "")


def test_gene_mask_is_case_insensitive(mg, genes_df):
    mask, success, message = mg.create_dataframe_gene_mask(genes_df, ["gene1", "GENE2"])
    assert genes_df.index[mask].tolist() == ["E1", "E2"]
    assert (success, message) == (1, "")


def test_gene_mask_picks_one_id_for_duplicated_symbols(mg, genes_df):
    mask, success, message = mg.create_dataframe_gene_mask(genes_df, ["Dup"])
    assert genes_df.index[mask].tolist() == ["E3"]
    assert success == 2 and "2 or more Ensembl IDs" in message and "Dup" in message


def test_gene_mask_reports_genes_not_found(mg, genes_df):
    mask, success, message = mg.create_dataframe_gene_mask(genes_df, ["Gene1", "Nope"])
    assert genes_df.index[mask].tolist() == ["E1"]
    assert success == 2 and "not found" in message and "Nope" in message


def test_gene_mask_errors_when_no_genes_found(mg, genes_df):
    with pytest.raises(mg.PlotError, match="None of the searched"):
        mg.create_dataframe_gene_mask(genes_df, ["Nope"])


def test_normalize_searched_genes(mg):
    normalized, found = mg.normalize_searched_genes(["Gene1", "GENE2", 5, "Other"], ["gene1", "gene2", "5", "missing"])
    assert normalized == ["Gene1", "GENE2", "5"]
    assert found == ["gene1", "gene2", "5"]


# ---------------------------------------------------------------- axis labels and titles

@pytest.fixture
def groups_df():
    return pd.DataFrame({"a": ["x", "x", "y", "x"], "b": ["p", "p", "q", "r"]})


def test_multicategory_labels_single_filter_is_flat(mg, groups_df):
    assert mg.create_multicategory_axis_labels(["a"], groups_df) == ["x", "x", "y", "x"]
    assert mg.create_multicategory_axis_labels(["a"], groups_df, return_unique=True) == ["x", "y"]


def test_multicategory_labels_two_filters_are_nested(mg, groups_df):
    assert mg.create_multicategory_axis_labels(["a", "b"], groups_df) == [["x", "x", "y", "x"], ["p", "p", "q", "r"]]
    assert mg.create_multicategory_axis_labels(["a", "b"], groups_df, return_unique=True) == [["x", "y", "x"], ["p", "q", "r"]]


def test_composite_index_column(mg):
    df = pd.DataFrame({"A": [1, 2], "B": ["x", "y"], "C": [3.5, None]})
    assert mg.create_composite_index_column(df, ["A", "B"]).tolist() == ["1;x", "2;y"]
    assert mg.create_composite_index_column(df, ["B", "C"]).tolist() == ["x;3.5", "y;nan"]


@pytest.mark.parametrize("filters, expected", [
    ([None], "Genes"),
    (["cell_type"], "Genes grouped by cell type"),
    (["cell_type", "time_point"], "Genes grouped by cell type and time point"),
])
def test_build_violin_x_title(mg, filters, expected):
    assert mg.build_violin_x_title(filters) == expected


def test_truncate_ticktext(mg):
    ticks, mapping = mg._truncate_ticktext(["short", "averyverylongname", "averyverylongother", 12345678901])
    assert ticks == ["short", "averyve...", "averyve...~1", "1234567..."]
    assert mapping == {"averyve...": "averyverylongname", "averyve...~1": "averyverylongother", "1234567...": "12345678901"}


@pytest.mark.parametrize("empty", [None, []])
def test_truncate_ticktext_empty(mg, empty):
    assert mg._truncate_ticktext(empty) == (None, {})


# ---------------------------------------------------------------- heatmap domains and palettes

def test_x_domains_with_dendrogram(mg):
    domains = mg.calculate_x_domains(show_dendrogram=True, gene_domain_ratio=0.5)
    assert domains["heatmap"] == pytest.approx([0.2, 0.9])
    assert domains["dendro"]["left"] == pytest.approx([0.02, 0.2])
    assert domains["dendro"]["top"] == pytest.approx([0.2, 0.55])


def test_x_domains_without_dendrogram(mg):
    domains = mg.calculate_x_domains(show_dendrogram=False)
    assert domains["heatmap"] == pytest.approx([0.05, 0.9])
    assert domains["dendro"]["top"] == pytest.approx(domains["heatmap"])


def test_y_domains_with_dendrogram(mg):
    domains = mg.calculate_y_domains(show_dendrogram=True, gene_domain_ratio=0.5)
    assert domains["heatmap"] == pytest.approx([0.1, 0.77])
    assert domains["dendro"]["left"] == pytest.approx([0.1, 0.435])
    assert domains["dendro"]["top"] == pytest.approx([0.77, 0.95])


def test_y_domains_without_dendrogram(mg):
    domains = mg.calculate_y_domains(show_dendrogram=False)
    assert domains["heatmap"] == pytest.approx([0.1, 0.95])
    assert domains["dendro"]["left"] == pytest.approx(domains["heatmap"])


def test_categorical_palette_cycles_through_four(mg):
    palettes = [mg.get_categorical_palette(i) for i in range(5)]
    assert palettes[0] == mg.cc.glasbey_dark and palettes[1] == mg.cc.glasbey_cool
    assert palettes[4] is palettes[0]
    assert len({id(p) for p in palettes[:4]}) == 4


def test_categorical_palette_colorblind(mg):
    assert all(mg.get_categorical_palette(i, colorblind=True) is mg.CATEGORICAL_PALETTE for i in range(3))


# ---------------------------------------------------------------- colorscale bounds (#1166)

@pytest.mark.parametrize("center, vmin, vmax, expected", [
    (False, None, None, (0, 5.0, None)),
    (False, 1.0, None, (1.0, 5.0, None)),
    (False, None, 3.0, (0, 3.0, None)),
    (True, None, None, (None, None, 0)),
    (True, None, 2.0, (-2.0, 2.0, None)),
    (True, -3.0, None, (-3.0, 3.0, None)),
    (True, 3.0, None, (-3.0, 3.0, None)),
    (True, -1.0, 4.0, (-1.0, 4.0, None)),
])
def test_heatmap_colorscale_bounds(mg, center, vmin, vmax, expected):
    assert mg._heatmap_colorscale_bounds(5.0, center, vmin, vmax) == expected


@pytest.mark.parametrize("center, vmin, vmax", [
    (False, 6.0, None),     # above the data maximum
    (False, None, -1.0),    # below the default minimum of 0
    (True, 1.0, 4.0),       # both on one side of zero
    (True, None, 0.0),
])
def test_heatmap_colorscale_bounds_errors(mg, center, vmin, vmax):
    with pytest.raises(mg.PlotError):
        mg._heatmap_colorscale_bounds(5.0, center, vmin, vmax)


def test_dot_plot_colorscale_bounds(mg):
    df = pd.DataFrame({
        "gene_symbol": ["G1", "G1", "G2", "G2"],
        "cluster": ["a", "b", "a", "b"],
        "mean": [1.0, 3.0, 7.0, 15.0],
        "count": [2, 2, 2, 2],
        "percent": [50.0, 100.0, 50.0, 100.0],
    })
    marker = mg.create_dot_plot(df, ["cluster"]).data[0].marker
    assert (marker.cmin, marker.cmax) == (0, 4.0)    # log2(15 + 1)

    marker = mg.create_dot_plot(df, ["cluster"], vmin=1.0, vmax=2.0).data[0].marker
    assert (marker.cmin, marker.cmax) == (1.0, 2.0)
