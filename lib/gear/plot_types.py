"""
plot_types.py - Lookups into the shared plot-type registry, www/plot_types.json.

The registry lists every display plot type with its family (which backend draws it), legacy
aliases, and layout hints. The browser reads the same file, so adding a plot type starts there
instead of in hard-coded lists scattered across the code.

The file is read on first use, not at import time: the upload-listener containers import geardb
(and so this module) but do not mount www/.
"""

import functools
import json
from pathlib import Path

REGISTRY_PATH = Path(__file__).resolve().parents[2] / "www" / "plot_types.json"

FAMILIES = {"plotly", "scanpy", "mg_plotly", "mg_scanpy", "svg", "gosling", "spatial"}
MULTIGENE_FAMILIES = {"mg_plotly", "mg_scanpy"}

# API path (appended to /api/plot/<dataset_id>) that draws a family's preview image.
# SVG previews are drawn from the uploaded SVG file instead, and gosling/spatial displays have no preview.
PREVIEW_ENDPOINTS = {
    "plotly": "",
    "mg_plotly": "/mg_plotly",
    "scanpy": "/tsne",
    "mg_scanpy": "/mg_tsne",
}


@functools.cache
def load_registry() -> dict[str, dict]:
    """Return the registry's plot types, keyed by canonical name."""
    try:
        with open(REGISTRY_PATH) as fh:
            return json.load(fh)["plot_types"]
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Plot-type registry not found at {REGISTRY_PATH}") from e


def _alias_map() -> dict[str, str]:
    """Map every canonical name and alias to its canonical name. Cheap to rebuild, so not cached separately."""
    lookup = {}
    for name, entry in load_registry().items():
        lookup[name] = name
        for alias in entry.get("aliases", []):
            lookup[alias] = name
    return lookup


def normalize(plot_type: str | None) -> str | None:
    """Return the canonical name for a plot type or one of its aliases, ignoring case and whitespace.

    Returns None for an empty or unknown plot type.
    """
    if not plot_type:
        return None
    return _alias_map().get(plot_type.strip().lower())


def _entry(plot_type: str | None) -> dict:
    name = normalize(plot_type)
    return load_registry()[name] if name else {}


def family(plot_type: str | None) -> str | None:
    """Return the family (plotly, scanpy, mg_plotly, ...) of a plot type, or None if unknown."""
    return _entry(plot_type).get("family")


def is_multigene(plot_type: str | None) -> bool:
    """Return True if the plot type shows several genes at once."""
    return family(plot_type) in MULTIGENE_FAMILIES


def is_full_width(plot_type: str | None) -> bool:
    """Return True if auto-arranged layouts should give this plot type a full-width tile."""
    return bool(_entry(plot_type).get("full_width", False))


def embedding_basis(plot_type: str | None) -> str | None:
    """Return the embedding (tsne, umap or pca) a scanpy-family plot type draws, or None."""
    return _entry(plot_type).get("basis")


def preview_endpoint(plot_type: str | None) -> str | None:
    """Return the API path suffix that draws this plot type's preview image, or None if the API does not draw it."""
    return PREVIEW_ENDPOINTS.get(family(plot_type))


def names(family_name: str | None = None) -> list[str]:
    """Return canonical plot type names, optionally only those of one family."""
    return [name for name, entry in load_registry().items() if family_name in (None, entry["family"])]
