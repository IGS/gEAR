"""
cgi_wb_real_analysis.py - The real gear.analysis, moved into a temporary tree.

root_dir becomes $GEAR_TEST_ANALYSIS_ROOT, so primary/public/user_saved analyses live under
<root>/www/..., and the user_unsaved area (/tmp/<session_id>/... in production) becomes
<root>/tmp/<session_id>/... Used by the stored-analysis CGI tests, which depend on the real
Analysis/AnalysisCollection behavior (pipeline JSON paths, from_json, vetting).
"""

import importlib.util
import os
import sys
from pathlib import Path

_REAL = Path(__file__).resolve().parents[3] / "lib" / "gear" / "analysis.py"
_spec = importlib.util.spec_from_file_location("_real_gear_analysis", _REAL)
_real = importlib.util.module_from_spec(_spec)
sys.modules["_real_gear_analysis"] = _real
_spec.loader.exec_module(_real)

_ROOT = Path(os.environ["GEAR_TEST_ANALYSIS_ROOT"])
_real.root_dir = _ROOT


def _out_of_tmp(path):
    """Move a /tmp/... path (the user_unsaved area) into <root>/tmp/...; paths already under root are kept."""
    path = Path(path)
    if path.is_relative_to(_ROOT) or not path.is_relative_to("/tmp"):
        return path
    return _ROOT / "tmp" / path.relative_to("/tmp")


_base_path = _real.Analysis.base_path.fget
_parent_path_by_type = _real.Analysis._parent_path_by_type
_real.Analysis.base_path = property(lambda self: _out_of_tmp(_base_path(self)))
_real.Analysis._parent_path_by_type = lambda self, atype=None: _out_of_tmp(_parent_path_by_type(self, atype))

globals().update({name: value for name, value in vars(_real).items() if not name.startswith("__")})
