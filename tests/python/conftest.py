"""
Shared pytest setup for the gEAR server-side test suite.

- Puts lib/, www/api/ and this directory on sys.path.
- Replaces gear.db.MySQLDB.connect with a fake connection before anything imports geardb:
  importing geardb opens a connection (LayoutCollection._cnx is created at class definition),
  and CI has neither MySQL nor a gear.ini.
- Provides fixtures for upload staging areas, example datasets and temporary gene list files.
"""

import shutil
import sys
import uuid
from pathlib import Path
from unittest import mock

import pytest

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parents[1]

for path in (REPO_ROOT / "lib", REPO_ROOT / "www" / "api", TESTS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import gear.db  # noqa: E402  (must be patched before geardb is imported anywhere)

gear.db.MySQLDB.connect = lambda self: mock.MagicMock(name="fake-mysql-connection")

WWW_DIR = REPO_ROOT / "www"
UPLOADS_DIR = WWW_DIR / "uploads" / "files"
CARTS_DIR = WWW_DIR / "carts"
USER_TEMPLATES = WWW_DIR / "user_templates"


@pytest.fixture
def example_mex_tar() -> Path:
    """The MEX example offered on the upload page (a plain, uncompressed tar)."""
    return USER_TEMPLATES / "example_mex.tar"


@pytest.fixture
def example_3tab_tar_gz() -> Path:
    """The 3-tab example offered on the upload page (a gzipped tar)."""
    return USER_TEMPLATES / "example_3tab.tar.gz"


@pytest.fixture
def upload_session():
    """
    A unique session directory under www/uploads/files (gitignored), removed afterwards.

    Yields (session_id, path). CGIs that locate uploads relative to their own file use this area.
    """
    session_id = f"pytest-{uuid.uuid4().hex[:12]}"
    path = UPLOADS_DIR / session_id
    path.mkdir(parents=True)
    try:
        yield session_id, path
    finally:
        shutil.rmtree(path, ignore_errors=True)


@pytest.fixture
def tmp_cart():
    """Write a weighted gene list file to www/carts (gitignored) and remove it afterwards."""
    created = []

    def _write(content: str) -> str:
        share_id = f"pytest{uuid.uuid4().hex[:10]}"
        path = CARTS_DIR / f"cart.{share_id}.tab"
        path.write_text(content)
        created.append(path)
        return share_id

    yield _write
    for path in created:
        path.unlink(missing_ok=True)
