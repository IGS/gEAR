"""Colorblind mode (issue #324): shared palettes, palette sampling, and the workbench/Gosling hooks."""

import importlib
import sys
import types
from pathlib import Path

import pytest

from gear import colorblind
from helpers.cgi_harness import run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"
WORKBENCH_ANALYSIS = {"gear.analysis": FAKES / "workbench_analysis.py"}
DB = {"sessions": {"owner": 1}, "datasets": [{"id": "DS1", "owner_id": 1}]}


# ---------------------------------------------------------------- lib/gear/colorblind.py

def test_categorical_colors_are_distinct_and_repeat_past_the_palette():
    palette_size = len(colorblind.CATEGORICAL_PALETTE)
    colors = colorblind.categorical_colors(palette_size + 2)
    assert len(set(colors[:palette_size])) == palette_size
    assert colors[palette_size:] == colors[:2]


@pytest.mark.parametrize("value, expected", [
    (True, True), ("true", True), ("1", True), ("True", True),
    (False, False), ("false", False), ("", False), (None, False), ("null", False),
])
def test_is_enabled(value, expected):
    assert colorblind.is_enabled(value) is expected


def test_remove_colorblind_copies_only_removes_named_copies(tmp_path):
    for name in ("tsne.png", "tsne_colorblind.png", "umap_colorblind.png"):
        (tmp_path / name).write_bytes(b"png")
    colorblind.remove_colorblind_copies(str(tmp_path), ["tsne", "missing"])
    assert sorted(p.name for p in tmp_path.iterdir()) == ["tsne.png", "umap_colorblind.png"]


# ---------------------------------------------------------------- mg_plotting palette sampling

@pytest.fixture
def mg_plotting(monkeypatch):
    """Import gear.mg_plotting, standing in for diffxpy (only used for differential expression) if it is missing."""
    pytest.importorskip("plotly")
    try:
        import diffxpy.api  # noqa: F401
    except ImportError:
        diffxpy = types.ModuleType("diffxpy")
        diffxpy.api = types.ModuleType("diffxpy.api")
        monkeypatch.setitem(sys.modules, "diffxpy", diffxpy)
        monkeypatch.setitem(sys.modules, "diffxpy.api", diffxpy.api)
    monkeypatch.delitem(sys.modules, "gear.mg_plotting", raising=False)
    return importlib.import_module("gear.mg_plotting")


def test_get_discrete_colors_samples_named_colorscales(mg_plotting):
    colors = mg_plotting.get_discrete_colors(["a", "b", "c"], "viridis")
    assert colors is not None and len(colors) == 3 and len(set(colors)) == 3
    assert mg_plotting.get_discrete_colors(["a", "b", "c"], "viridis", reverse_colorscale=True) == colors[::-1]



def test_colorblind_swatch_gives_palette_colors_in_order(mg_plotting):
    assert mg_plotting.get_discrete_colors(["a", "b", "c"], mg_plotting.COLORBLIND_SWATCH)[:3] == colorblind.categorical_colors(3)


def test_quadrant_plot_colorblind_swatch(mg_plotting):
    import pandas as pd
    df = pd.DataFrame({"s1_c_log2FC": [1.0, -1.0], "s2_c_log2FC": [1.0, -1.0],
                       "ensm_id": ["E1", "E2"], "gene_symbol": ["G1", "G2"]})
    fig = mg_plotting.create_quadrant_plot(df, "ctrl", "c1", "c2", colorscale=mg_plotting.COLORBLIND_SWATCH)
    # UP/UP and DOWN/DOWN take the first two palette colors
    assert [trace.marker.color for trace in fig.data] == colorblind.categorical_colors(2)

def test_quadrant_plot_uses_sampled_colorscale(mg_plotting):
    import pandas as pd
    df = pd.DataFrame({"s1_c_log2FC": [1.0, -1.0], "s2_c_log2FC": [1.0, -1.0],
                       "ensm_id": ["E1", "E2"], "gene_symbol": ["G1", "G2"]})
    fig = mg_plotting.create_quadrant_plot(df, "ctrl", "c1", "c2", colorscale="viridis")
    expected = mg_plotting.px.colors.sample_colorscale(mg_plotting.px.colors.get_colorscale("viridis"), 8)[:2]
    assert [trace.marker.color for trace in fig.data] == expected

# ---------------------------------------------------------------- Gosling spec

@pytest.fixture(scope="module")
def gosling_spec():
    pytest.importorskip("gosling")
    from resources import gosling_spec
    return gosling_spec


def _strand_colors(view):
    return [track["color"]["range"] for track in view.to_dict()["tracks"] if track.get("color", {}).get("field") == "strand"]


@pytest.mark.parametrize("colorblind_mode", [False, True])
def test_annotation_strand_colors(gosling_spec, colorblind_mode):
    view = gosling_spec.build_bed_annotation_tracks("mm10", colorblind=colorblind_mode)
    expected = colorblind.categorical_colors(2) if colorblind_mode else ["darkblue", "darkred"]
    ranges = _strand_colors(view)
    assert ranges and all(r == expected for r in ranges)


@pytest.mark.parametrize("colorblind_mode, colorscale", [(False, "bupu"), (True, "cividis")])
def test_hic_colorscale(gosling_spec, colorblind_mode, colorscale):
    spec = gosling_spec.HiCSpec(data_url="https://example.org/a.mcool", colorblind=colorblind_mode)
    track = spec.get_encoding(100, 100).to_dict()
    assert track["color"]["range"] == colorscale



def _track_colors(view):
    """Single colors of every data track in a rendered panel, in order (multiWig members included)."""
    def walk(node):
        if isinstance(node, dict):
            color = node.get("color")
            if isinstance(color, dict) and "value" in color:
                yield color["value"]
            for value in node.values():
                yield from walk(value)
        elif isinstance(node, list):
            for item in node:
                yield from walk(item)
    return list(walk(view.to_dict()))


HUB_TRACKS = [
    {"type": "bigWig", "track": "a", "bigDataUrl": "https://example.org/a.bw", "color": "rgb(255,0,0)"},
    {"type": "bigWig", "track": "b", "bigDataUrl": "https://example.org/b.bw", "color": "rgb(0,255,0)"},
    {"track": "overlay", "container": "multiWig"},
    {"type": "bigWig", "track": "c", "parent": "overlay", "bigDataUrl": "https://example.org/c.bw", "color": "rgb(255,0,0)"},
    {"type": "bigWig", "track": "d", "parent": "overlay", "bigDataUrl": "https://example.org/d.bw", "color": "rgb(0,255,0)"},
]


def test_hub_track_colors_kept_in_normal_mode(gosling_spec):
    view, _ = gosling_spec.build_gosling_tracks([], HUB_TRACKS)
    assert set(_track_colors(view)) == {"rgb(255,0,0)", "rgb(0,255,0)"}


def test_hub_track_colors_replaced_in_colorblind_mode(gosling_spec):
    view, _ = gosling_spec.build_gosling_tracks([], HUB_TRACKS, colorblind=True)
    colors = _track_colors(view)
    # Every data track gets its own palette color, skipping the two used by the annotation strands
    assert len(colors) == 4 and len(set(colors)) == 4
    assert set(colors) <= set(colorblind.categorical_colors(len(HUB_TRACKS) + 2)[2:])


# ---------------------------------------------------------------- heatmap clusterbars

@pytest.mark.parametrize("colorblind_mode", [False, True])
def test_clusterbar_palette(mg_plotting, colorblind_mode):
    import pandas as pd
    import plotly.graph_objects as go
    obs_columns = pd.MultiIndex.from_tuples([("a", "x"), ("b", "y")], names=["cluster", "batch"])
    fig = go.Figure()
    mg_plotting.add_clusterbars(fig, obs_columns, ["cluster", "batch"], 1.0, colorblind=colorblind_mode)
    first_colors = [trace.colorscale[0][1] for trace in fig.data]
    if colorblind_mode:
        assert first_colors == [colorblind.CATEGORICAL_PALETTE[0]] * 2
    else:
        assert first_colors == [mg_plotting.cc.glasbey_dark[0], mg_plotting.cc.glasbey_cool[0]]

# ---------------------------------------------------------------- workbench images

class TestAnalysisImage:
    """get_analysis_image.cgi serves a step's colorblind copy only to colorblind viewers."""

    def run(self, tmp_path, **query):
        query = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "user_unsaved",
                 "session_id": "owner", "analysis_name": "tsne", **query}
        return run_cgi("get_analysis_image.cgi", db=DB, query=query, fake_modules=WORKBENCH_ANALYSIS,
                       env={"GEAR_TEST_ANALYSIS_DIR": str(tmp_path)})

    @pytest.fixture
    def figures(self, tmp_path):
        figures = tmp_path / "user_unsaved" / "figures"
        figures.mkdir(parents=True)
        (figures / "tsne.png").write_bytes(b"normal")
        (figures / "tsne_colorblind.png").write_bytes(b"colorblind")
        return figures

    @pytest.mark.parametrize("flag, expected", [("true", b"colorblind"), ("false", b"normal"), (None, b"normal")])
    def test_serves_copy_only_in_colorblind_mode(self, tmp_path, figures, flag, expected):
        query = {} if flag is None else {"colorblind_mode": flag}
        assert self.run(tmp_path, **query).body.strip() == expected

    def test_falls_back_when_step_has_no_copy(self, tmp_path, figures):
        (figures / "pca.png").write_bytes(b"pca")
        assert self.run(tmp_path, analysis_name="pca", colorblind_mode="true").body.strip() == b"pca"

    def test_rejects_path_traversal(self, tmp_path, figures):
        result = self.run(tmp_path, analysis_name="../../../../etc/passwd", colorblind_mode="true")
        assert b"root:" not in result.body



class TestStepCopies:
    """Plotting steps save a colorblind copy only in colorblind mode, and remove a stale one otherwise."""

    TSNE_QUERY = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "user_unsaved", "session_id": "owner",
                  "n_pcs": "2", "n_neighbors": "5", "random_state": "0", "genes_to_color": "G1", "use_scaled": "false",
                  "compute_neighbors": "0", "compute_tsne": "0", "compute_umap": "0", "plot_tsne": "1", "plot_umap": "0"}
    CLUSTER_QUERY = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "user_unsaved", "session_id": "owner",
                     "resolution": "1.0", "compute_clusters": "false", "cluster_info": "[]",
                     "plot_tsne": "1", "plot_umap": "0"}

    @pytest.fixture(autouse=True)
    def needs_scanpy(self):
        pytest.importorskip("scanpy")

    def run(self, tmp_path, script, query):
        (tmp_path / "user_unsaved").mkdir(exist_ok=True)
        result = run_cgi(script, db=DB, query=query, fake_modules=WORKBENCH_ANALYSIS,
                         env={"GEAR_TEST_ANALYSIS_DIR": str(tmp_path)})
        assert result.returncode == 0, result.stderr[-3000:]
        assert result.json()["success"] == 1, result.stderr[-3000:]
        return sorted(p.name for p in (tmp_path / "user_unsaved" / "figures").iterdir())

    @pytest.mark.parametrize("script, query, name", [
        ("h5ad_generate_tsne.cgi", TSNE_QUERY, "tsne"),
        ("h5ad_generate_clusters.cgi", CLUSTER_QUERY, "tsne_clustering"),
    ])
    def test_copy_made_then_removed(self, tmp_path, script, query, name):
        assert self.run(tmp_path, script, {**query, "colorblind_mode": "true"}) == \
            [f"{name}.png", f"{name}_colorblind.png"]
        assert self.run(tmp_path, script, {**query, "colorblind_mode": "false"}) == [f"{name}.png"]

    def test_marker_gene_visualization_copies(self, tmp_path):
        query = {"dataset_id": "DS1", "analysis_id": "A1", "analysis_type": "user_unsaved", "session_id": "owner",
                 "marker_genes": '["G1", "G2"]'}
        script = "h5ad_generate_marker_gene_visualization.cgi"
        assert self.run(tmp_path, script, {**query, "colorblind_mode": "true"}) == [
            "dotplot_goi.png", "dotplot_goi_colorblind.png", "stacked_violin_goi.png", "stacked_violin_goi_colorblind.png"]
        assert self.run(tmp_path, script, {**query, "colorblind_mode": "false"}) == [
            "dotplot_goi.png", "stacked_violin_goi.png"]

    def test_primary_filter_copies(self, tmp_path):
        query = {"dataset_id": "DS1", "analysis_type": "primary", "session_id": "owner",
                 "filter_cells_lt_n_genes": "", "filter_cells_gt_n_genes": "",
                 "filter_genes_lt_n_cells": "", "filter_genes_gt_n_cells": ""}
        script = "h5ad_apply_primary_filter.cgi"
        assert self.run(tmp_path, script, {**query, "colorblind_mode": "true"}) == [
            "highest_expr_genes.png", "highest_expr_genes_colorblind.png"]
        assert self.run(tmp_path, script, {**query, "colorblind_mode": "false"}) == ["highest_expr_genes.png"]

    def test_tsne_without_genes_makes_no_copy(self, tmp_path):
        query = {**self.TSNE_QUERY, "genes_to_color": "", "colorblind_mode": "true"}
        assert self.run(tmp_path, "h5ad_generate_tsne.cgi", query) == ["tsne.png"]

# ---------------------------------------------------------------- h5ad_find_marker_genes.cgi

BICCN_DATASET = "4fbd43e2-f301-42d8-a61d-a5dd5bf720e7"


def test_marker_genes_uses_biccn_cluster_column_for_primary_analysis(tmp_path):
    pytest.importorskip("scanpy")
    db = {"sessions": {"owner": 1}, "datasets": [{"id": BICCN_DATASET, "owner_id": 1}]}
    result = run_cgi("h5ad_find_marker_genes.cgi", db=db, fake_modules=WORKBENCH_ANALYSIS,
                     env={"GEAR_TEST_ANALYSIS_DIR": str(tmp_path)},
                     query={"dataset_id": BICCN_DATASET, "analysis_type": "primary", "session_id": "owner",
                            "n_genes": "3", "compute_marker_genes": "true"})
    assert result.returncode == 0, result.stderr[-3000:]
    assert "joint_cluster_round4_annot" in str(result.json()), result.stderr[-3000:]
    # Results are still written to the user's unsaved analysis, not the primary one
    assert (tmp_path / "user_unsaved" / f"{BICCN_DATASET}.h5ad").exists()
    assert not (tmp_path / "primary").exists()
