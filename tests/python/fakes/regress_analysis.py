"""
The real gear.analysis with its root_dir moved to $REGRESS_ROOT_DIR, so Analysis.base_path
(primary -> <root>/www/datasets, user_saved -> <root>/www/analyses/by_user/...) points into a
temporary tree. Used by delete_dataset_analysis.cgi regression tests.
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

_real.root_dir = Path(os.environ["REGRESS_ROOT_DIR"])

Analysis = _real.Analysis
