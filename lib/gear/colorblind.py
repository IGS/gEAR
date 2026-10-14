"""
colorblind.py - Palettes used when a user has colorblind mode turned on (issue #324).

Colorblind mode only changes what the current viewer sees; callers must not save these
palettes into displays or images that other users load.
"""

import colorcet as cc

# Continuous scales (e.g. gene expression); matches tsne_data.py and plotly_data.py
CONTINUOUS_CMAP = "cividis_r"

# Many distinct colors from the blue/green/purple range, avoiding red-green pairs
CATEGORICAL_PALETTE = cc.glasbey_cool


def categorical_colors(n):
    """Return n hex colors for categorical data, repeating the palette if n exceeds its length."""
    return [CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)] for i in range(n)]


def variant_name(name):
    """Return the name of the colorblind copy of a saved analysis figure (e.g. 'tsne' -> 'tsne_colorblind')."""
    return f"{name}_colorblind"


def is_enabled(value):
    """Interpret a colorblind_mode request value (form fields arrive as strings like 'true' or '1')."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("1", "true", "yes")


def remove_colorblind_copies(figures_dir, names):
    """
    Delete stale colorblind copies of the given figures (e.g. after a step is re-run without colorblind
    mode), so a colorblind viewer never sees an image left over from an earlier run.
    """
    import os

    for name in names:
        path = os.path.join(figures_dir, f"{variant_name(name)}.png")
        if os.path.exists(path):
            os.remove(path)
