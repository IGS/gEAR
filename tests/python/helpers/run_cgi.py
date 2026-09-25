"""
run_cgi.py - Run a gEAR CGI script as __main__ with test doubles installed first.

Usage (normally via cgi_harness.run_cgi): python run_cgi.py <path/to/script.cgi>

Before the script runs, the fake geardb (helpers/fake_geardb.py) is placed in sys.modules, so it
is used even by scripts that put the real lib/ first on sys.path. Extra modules can be replaced
with GEAR_FAKE_MODULES, a JSON object of {module name: path to a .py file}.
"""

import importlib.util
import json
import os
import runpy
import sys
from pathlib import Path

HELPERS_DIR = Path(__file__).resolve().parent
REPO_ROOT = HELPERS_DIR.parents[2]


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    script = Path(sys.argv[1]).resolve()

    # The real gear package stays importable (for modules that aren't faked)
    sys.path.append(str(REPO_ROOT / "lib"))

    _load_module("geardb", HELPERS_DIR / "fake_geardb.py")
    for name, path in json.loads(os.environ.get("GEAR_FAKE_MODULES", "{}")).items():
        _load_module(name, path)

    sys.argv = [str(script)]
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()
