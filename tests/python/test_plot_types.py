"""
The shared plot-type registry (www/plot_types.json) and its Python helper, gear.plot_types (#1380).

The expected values below are copied from the hard-coded lists the registry replaces, so these
tests show the move did not change behavior. The sync tests keep the registry, the API's
available display types and the curator plot-type menus from drifting apart.
"""

import json
import re
import types
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest
from flask import Flask
from flask_restful import Api
from gear import plot_types
from resources import available_display_types

WWW = Path(__file__).resolve().parents[2] / "www"

SINGLE_GENE_TYPES = {
    "bar", "line", "scatter", "scatter3d", "violin", "tsne/umap_dynamic",
    "tsne_static", "umap_static", "pca_static", "svg", "gosling", "spatial_panel",
}
# Was multi_plot_types in geardb.LayoutDisplay.get_is_multigene
MULTIGENE_TYPES = {
    "dotplot", "heatmap", "mg_violin", "quadrant", "volcano", "mg_tsne_static", "mg_umap_static", "mg_pca_static",
}
# Was geardb.EMBEDDING_PLOT_TYPES, minus scatter3d: its boxy 3D plots read fine in a narrower tile
FULL_WIDTH_TYPES = {
    "tsne_static", "umap_static", "pca_static", "tsne/umap_dynamic",
    "mg_tsne_static", "mg_umap_static", "mg_pca_static",
}
LEGACY_ALIASES = {"tsne": "tsne_static", "tsne_dynamic": "tsne/umap_dynamic", "tsne_dyna": "tsne/umap_dynamic"}
# Was tsne_data.PLOT_TYPE_TO_BASIS
EMBEDDING_BASIS = {
    "tsne_static": "tsne", "tsne": "tsne", "umap_static": "umap", "pca_static": "pca",
    "mg_tsne_static": "tsne", "mg_umap_static": "umap", "mg_pca_static": "pca",
}


@pytest.fixture
def registry():
    return json.loads((WWW / "plot_types.json").read_text())["plot_types"]


# ---------------------------------------------------------------- registry file

class TestRegistrySchema:
    def test_lists_every_plot_type(self, registry):
        assert set(registry) == SINGLE_GENE_TYPES | MULTIGENE_TYPES

    def test_entries_have_a_label_and_known_family(self, registry):
        for name, entry in registry.items():
            assert entry["family"] in plot_types.FAMILIES, name
            assert isinstance(entry["label"], str) and entry["label"], name

    def test_aliases_are_unique_and_never_a_plot_type_name(self, registry):
        aliases = [alias for entry in registry.values() for alias in entry.get("aliases", [])]
        assert len(aliases) == len(set(aliases))
        assert not set(aliases) & set(registry)

    def test_only_embedding_families_have_a_basis(self, registry):
        for name, entry in registry.items():
            assert ("basis" in entry) == (entry["family"] in {"scanpy", "mg_scanpy"}), name


# ---------------------------------------------------------------- helper lookups

class TestLookups:
    @pytest.mark.parametrize("name", sorted(SINGLE_GENE_TYPES | MULTIGENE_TYPES))
    def test_canonical_names_normalize_to_themselves(self, name):
        assert plot_types.normalize(name) == name

    @pytest.mark.parametrize("alias, name", LEGACY_ALIASES.items())
    def test_legacy_aliases(self, alias, name):
        assert plot_types.normalize(alias) == name

    @pytest.mark.parametrize("value", ["  Bar ", "TSNE_STATIC", "Tsne"])
    def test_normalize_ignores_case_and_whitespace(self, value):
        assert plot_types.normalize(value) in {"bar", "tsne_static"}

    @pytest.mark.parametrize("value", [None, "", "pie", "epiviz", "contour"])
    def test_unknown_plot_types(self, value):
        assert plot_types.normalize(value) is None
        assert plot_types.family(value) is None
        assert plot_types.is_multigene(value) is False
        assert plot_types.is_full_width(value) is False
        assert plot_types.embedding_basis(value) is None
        assert plot_types.preview_endpoint(value) is None

    @pytest.mark.parametrize("name", sorted(SINGLE_GENE_TYPES | MULTIGENE_TYPES | set(LEGACY_ALIASES)))
    def test_is_multigene(self, name):
        assert plot_types.is_multigene(name) is (name in MULTIGENE_TYPES)

    @pytest.mark.parametrize("name", sorted(SINGLE_GENE_TYPES | MULTIGENE_TYPES | set(LEGACY_ALIASES)))
    def test_is_full_width(self, name):
        """Legacy names now get the same width as their current names (they were missed before)."""
        assert plot_types.is_full_width(name) is (plot_types.normalize(name) in FULL_WIDTH_TYPES)

    @pytest.mark.parametrize("name", sorted(SINGLE_GENE_TYPES | MULTIGENE_TYPES | set(LEGACY_ALIASES)))
    def test_embedding_basis(self, name):
        assert plot_types.embedding_basis(name) == EMBEDDING_BASIS.get(name, EMBEDDING_BASIS.get(plot_types.normalize(name)))

    @pytest.mark.parametrize("name, endpoint", [
        ("bar", ""), ("scatter3d", ""), ("tsne_dynamic", ""), ("violin", ""),
        ("dotplot", "/mg_plotly"), ("quadrant", "/mg_plotly"),
        ("tsne", "/tsne"), ("pca_static", "/tsne"),
        ("mg_pca_static", "/mg_tsne"),
        ("svg", None), ("gosling", None), ("spatial_panel", None),
    ])
    def test_preview_endpoint(self, name, endpoint):
        assert plot_types.preview_endpoint(name) == endpoint

    def test_names_by_family(self):
        assert set(plot_types.names("mg_scanpy")) == {"mg_tsne_static", "mg_umap_static", "mg_pca_static"}
        assert set(plot_types.names()) == SINGLE_GENE_TYPES | MULTIGENE_TYPES

    def test_missing_registry_file_gives_a_clear_error(self, monkeypatch, tmp_path):
        monkeypatch.setattr(plot_types, "REGISTRY_PATH", tmp_path / "plot_types.json")
        plot_types.load_registry.cache_clear()
        try:
            with pytest.raises(FileNotFoundError, match="plot_types.json"):
                plot_types.normalize("bar")
        finally:
            plot_types.load_registry.cache_clear()


# ---------------------------------------------------------------- sync with the API and the curator menus

def available_types(monkeypatch, tmp_path, resource):
    """Call an available-display-types resource for a dataset that supports every plot type."""
    n = 6
    obs = pd.DataFrame(
        {
            "cluster": pd.Categorical(["a", "b", "c"] * 2),
            "batch": pd.Categorical(["x", "y"] * 3),
            "sex": pd.Categorical(["f", "m", "f", "m", "m", "f"]),
            "time_point": np.arange(n, dtype=float),
            "score": np.linspace(0, 1, n),
        },
        index=[f"cell{i}" for i in range(n)],
    )
    adata = anndata.AnnData(X=np.ones((n, 2), dtype=np.float32), obs=obs)
    adata.obsm["X_umap"] = np.zeros((n, 2))
    ds = types.SimpleNamespace(
        id="DS1", dtype="single-cell-rnaseq", get_file_path=lambda: str(tmp_path / "datasets" / "DS1.h5ad")
    )
    monkeypatch.setattr(available_display_types.geardb, "get_dataset_by_id", lambda *a, **k: ds)
    monkeypatch.setattr(available_display_types, "get_adata_shadow", lambda *a, **k: adata)

    app = Flask(__name__)
    Api(app).add_resource(resource, "/types/<dataset_id>")
    resp = app.test_client().post("/types/DS1", json={"dataset_id": "DS1"})
    assert resp.status_code == 200, resp.get_data(as_text=True)
    return resp.get_json()


@pytest.mark.parametrize("resource, multigene", [
    (available_display_types.AvailableDisplayTypes, False),
    (available_display_types.MGAvailableDisplayTypes, True),
])
def test_available_display_types_are_registry_types(monkeypatch, tmp_path, resource, multigene):
    """Every plot type the API offers is a registry type of the right kind (single-gene or multigene)."""
    for name in available_types(monkeypatch, tmp_path, resource):
        assert plot_types.normalize(name) == name, name
        assert plot_types.is_multigene(name) is multigene, name


def curator_plot_type_options(page):
    """The plot-type menu of a curator page, as {option value: option text}, minus the "Select plot type" placeholder."""
    html = (WWW / page).read_text()
    menu = re.search(r'<select id="plot-type-select".*?</select>', html, re.DOTALL).group(0)
    options = dict(re.findall(r'<option value="([^"]+)"[^>]*>([^<]+)</option>', menu))
    options.pop("nope", None)
    return options


@pytest.mark.parametrize("page, families", [
    ("dataset_curator.html", {"plotly", "scanpy", "svg"}),
    ("multigene_curator.html", {"mg_plotly", "mg_scanpy"}),
])
def test_curator_menu_matches_registry(page, families):
    """Each curator offers exactly the registry types it can build. Gosling and spatial displays are made elsewhere."""
    offered = {plot_types.normalize(value) for value in curator_plot_type_options(page)}
    expected = {name for family in families for name in plot_types.names(family)}
    assert offered == expected


@pytest.mark.parametrize("page", ["dataset_curator.html", "multigene_curator.html"])
def test_curator_menu_text_matches_registry_labels(page, registry):
    """Labels shown elsewhere (explorer tags, display cards) match what curators pick from the menu."""
    for value, text in curator_plot_type_options(page).items():
        assert registry[plot_types.normalize(value)]["label"] == text.strip(), value
