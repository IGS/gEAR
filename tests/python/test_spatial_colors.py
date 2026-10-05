"""services/spatial/common.get_color_maps: per-viewer colorblind palettes for the spatial Panel viewer (issue #324)."""

import importlib.util
from pathlib import Path

import pytest

from gear import colorblind

pytestmark = pytest.mark.spatial

for module in ("datashader", "holoviews", "bokeh", "param", "colorcet", "matplotlib"):
    pytest.importorskip(module)

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
