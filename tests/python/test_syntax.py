"""Every server-side Python file (including CGI scripts) must at least parse."""

import ast
import warnings
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]


def _python_files():
    for folder in ("lib", "www/api", "listeners", "services"):
        yield from (p for p in (REPO_ROOT / folder).rglob("*.py") if "__pycache__" not in p.parts)
    for path in sorted((REPO_ROOT / "www" / "cgi").iterdir()):
        if path.suffix == ".py":
            yield path
        elif path.suffix == ".cgi" and "python" in path.read_text(errors="ignore").split("\n", 1)[0]:
            yield path


PYTHON_FILES = sorted(set(_python_files()))


@pytest.mark.parametrize("path", PYTHON_FILES, ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_parses(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", SyntaxWarning)   # e.g. invalid escapes in older regex strings
        ast.parse(path.read_text(), filename=str(path))
