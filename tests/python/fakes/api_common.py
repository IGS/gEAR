"""
api_common.py - Stand-in for www/api/resources/common.py, used only when the real module can't
be imported (it needs plotly and anndata "shadows", which the test requirements don't install).

Tests monkeypatch these functions on the resource module that imports them, so the bodies here
only fail loudly if a test forgets to.
"""


def _not_patched(*args, **kwargs):
    raise AssertionError("resources.common stub called; monkeypatch this function in the test")


get_adata_from_analysis = _not_patched
get_adata_shadow = _not_patched
get_adata_shadow_from_analysis = _not_patched
get_spatial_adata = _not_patched
