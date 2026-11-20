"""services/spatial/common.get_color_maps: per-viewer colorblind palettes for the spatial Panel viewer (issue #324)."""

import importlib.util
from pathlib import Path

import pytest

from gear import colorblind

pytestmark = pytest.mark.spatial

for module in ("datashader", "holoviews", "bokeh", "param", "colorcet", "matplotlib"):
    pytest.importorskip(module)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SPATIAL_COMMON = Path(__file__).resolve().parents[2] / "services" / "spatial" / "common.py"


@pytest.fixture(scope="module")
def spatial_common():
    spec = importlib.util.spec_from_file_location("spatial_common", SPATIAL_COMMON)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def df():
    # Colors as stored in the cached gene CSV
    return pd.DataFrame({"clusters": ["10", "2", "2", "1"], "colors": ["#aa0000", "#00aa00", "#00aa00", "#0000aa"]})


def test_normal_mode_uses_cached_colors(spatial_common, df):
    expression_cmap, cluster_cmap = spatial_common.get_color_maps(df, colorblind_mode=False)
    assert expression_cmap is spatial_common.cc.m_CET_L4_r
    assert cluster_cmap == {"10": "#aa0000", "2": "#00aa00", "1": "#0000aa"}


def test_colorblind_mode_uses_cividis_and_glasbey_cool(spatial_common, df):
    expression_cmap, cluster_cmap = spatial_common.get_color_maps(df, colorblind_mode=True)
    assert expression_cmap.name == colorblind.CONTINUOUS_CMAP
    # Clusters sorted numerically, colored in palette order, matching lib/gear/colorblind.py
    assert cluster_cmap == dict(zip(["1", "2", "10"], colorblind.categorical_colors(3)))


@pytest.mark.parametrize("cbar_min, cbar_max, expected", [
    (None, None, None),     # HoloViews picks the limits
    (None, 5.0, (0, 5.0)),  # the 98th-percentile default
    (1.0, 5.0, (1.0, 5.0)),
])
def test_expression_clim(spatial_common, cbar_min, cbar_max, expected):
    """#1166: user-defined bounds for the expression colorbar."""
    assert spatial_common.expression_clim(cbar_min, cbar_max) == expected


def test_expression_clim_without_upper_limit_leaves_it_to_the_data(spatial_common):
    low, high = spatial_common.expression_clim(2.0, None)
    assert low == 2.0 and np.isnan(high)
