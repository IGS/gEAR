"""services/spatial/common.apply_range: view ranges and flipped axes for the spatial Panel viewer (issue #1301)."""

import importlib.util
from pathlib import Path

import pytest

pytestmark = pytest.mark.spatial

for module in ("datashader", "holoviews", "bokeh", "param", "colorcet", "matplotlib"):
    pytest.importorskip(module)

from bokeh.models import Range1d  # noqa: E402

SPATIAL_COMMON = Path(__file__).resolve().parents[2] / "services" / "spatial" / "common.py"


@pytest.fixture(scope="module")
def spatial_common():
    spec = importlib.util.spec_from_file_location("spatial_common", SPATIAL_COMMON)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def full_range():
    # As HoloViews renders it: the full extent, with matching reset values
    return Range1d(start=0, end=100, reset_start=0, reset_end=100)


def test_no_bounds_keeps_range_and_reset(spatial_common):
    rng = full_range()
    spatial_common.apply_range(rng)
    assert (rng.start, rng.end) == (0, 100)
    assert (rng.reset_start, rng.reset_end) == (0, 100)


def test_no_bounds_fills_missing_reset_from_start_end(spatial_common):
    rng = Range1d(start=0, end=100)
    spatial_common.apply_range(rng)
    assert (rng.reset_start, rng.reset_end) == (0, 100)


def test_zoomed_bounds_keep_full_extent_for_reset(spatial_common):
    rng = full_range()
    spatial_common.apply_range(rng, (10, 50))
    assert (rng.start, rng.end) == (10, 50)
    assert (rng.reset_start, rng.reset_end) == (0, 100)


def test_flipped_bounds_flip_reset_too(spatial_common):
    """A flipped range carried over from the previous gene stays flipped after Reset."""
    rng = full_range()
    spatial_common.apply_range(rng, (50, 10))
    assert (rng.start, rng.end) == (50, 10)
    assert (rng.reset_start, rng.reset_end) == (100, 0)


def test_unflipped_bounds_unflip_a_flipped_reset(spatial_common):
    rng = Range1d(start=100, end=0, reset_start=100, reset_end=0)
    spatial_common.apply_range(rng, (10, 50))
    assert (rng.reset_start, rng.reset_end) == (0, 100)


SPATIAL_DIR = SPATIAL_COMMON.parent


@pytest.fixture(scope="module")
def panel_common(monkeypatch_module):
    pytest.importorskip("panel")
    pytest.importorskip("hvplot")
    # panel_common imports its sibling module as "common"
    monkeypatch_module.syspath_prepend(str(SPATIAL_DIR))
    import panel_common as module
    return module


@pytest.fixture(scope="module")
def monkeypatch_module():
    with pytest.MonkeyPatch.context() as mp:
        yield mp


class TestRestoredRangeReachesFrontendState:
    def test_restored_range_is_in_state_before_any_pan_or_zoom(self, panel_common):
        """Regression for #1301: a session that restores a carried-over view range (gene switch, expanded view)
        includes that range in the state sent to the browser right away. The HTML download and "save display"
        read that state, so without it they lost the pan/zoom/flip until the user panned again."""
        import json
        import types

        viewer = types.SimpleNamespace(
            settings=panel_common.Settings(dataset_id="DS1", filename="Actb.csv"),
            restore_xlim=(80.0, 20.0),  # flipped X
            restore_ylim=(10.0, 90.0),
            state_sync=types.SimpleNamespace(value=""),
        )
        figure = types.SimpleNamespace(x_range=Range1d(start=0, end=100), y_range=Range1d(start=0, end=100))

        panel_common.BaseSpatialViewer._apply_restored_range(viewer, figure)
        panel_common.BaseSpatialViewer._sync_state_to_js(viewer)

        state = json.loads(viewer.state_sync.value)
        range_keys = ("x_range_start", "x_range_end", "y_range_start", "y_range_end")
        assert {k: state.get(k) for k in range_keys} == {
            "x_range_start": 80.0, "x_range_end": 20.0, "y_range_start": 10.0, "y_range_end": 90.0,
        }
