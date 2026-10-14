"""lib/geardb.py helpers that don't need a database."""

import importlib

import pytest


def test_geardb_imports_without_mysql():
    """conftest.py replaces MySQLDB.connect, so importing geardb needs neither MySQL nor gear.ini."""
    geardb = importlib.import_module("geardb")
    assert hasattr(geardb, "Dataset")
    assert callable(geardb.get_dataset_by_id)


def test_importing_geardb_does_not_connect():
    """Only code that queries should open a connection (checked in a fresh interpreter)."""
    import subprocess
    import sys
    from pathlib import Path

    lib = Path(__file__).resolve().parents[2] / "lib"
    code = (
        "import sys; sys.path.insert(0, %r)\n"
        "import gear.db\n"
        "def refuse(self): raise AssertionError('connected at import')\n"
        "gear.db.MySQLDB.connect = refuse\n"
        "import geardb\n"
    ) % str(lib)
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, proc.stderr[-2000:]


class TestLayoutCollectionConnection:
    @pytest.fixture
    def connections(self, monkeypatch):
        import geardb
        from unittest import mock

        opened = []

        def connect(self):
            cnx = mock.MagicMock(name=f"cnx{len(opened)}")
            cnx.cursor.return_value.fetchall.return_value = []
            cnx.is_connected.return_value = True
            opened.append(cnx)
            return cnx

        monkeypatch.setattr(geardb.gear.db.MySQLDB, "connect", connect)
        monkeypatch.setattr(geardb.LayoutCollection, "_shared_cnx", None)
        return opened

    def test_collections_share_one_connection(self, connections):
        import geardb

        geardb.LayoutCollection()
        geardb.LayoutCollection().get_by_share_id("nope")
        assert len(connections) == 1
        assert connections[0].cursor.call_count >= 2

    def test_reconnects_after_the_server_drops_it(self, connections):
        import geardb

        geardb.LayoutCollection()
        connections[0].is_connected.return_value = False
        geardb.LayoutCollection()
        assert len(connections) == 2


class TestGetSourceArchivePath:
    @pytest.fixture
    def dataset(self, monkeypatch, tmp_path):
        import geardb

        ds = geardb.Dataset(id="ID", dtype="single-cell-rnaseq")
        monkeypatch.setattr(ds, "get_tarball_path", lambda: str(tmp_path / "ID.tar.gz"))
        return ds

    def test_no_archive(self, dataset):
        assert dataset.get_source_archive_path() is None

    @pytest.mark.parametrize("name", ["ID.tar.gz", "ID.tar", "ID.zip"])
    def test_archive_with_each_extension(self, dataset, tmp_path, name):
        (tmp_path / name).write_bytes(b"archive")
        assert dataset.get_source_archive_path() == str(tmp_path / name)

    def test_tar_gz_preferred_over_other_extensions(self, dataset, tmp_path):
        for name in ("ID.zip", "ID.tar", "ID.tar.gz"):
            (tmp_path / name).write_bytes(b"archive")
        assert dataset.get_source_archive_path() == str(tmp_path / "ID.tar.gz")

    def test_other_dataset_archives_are_ignored(self, dataset, tmp_path):
        (tmp_path / "OTHER.tar.gz").write_bytes(b"archive")
        (tmp_path / "ID.h5ad").write_bytes(b"not an archive")
        assert dataset.get_source_archive_path() is None

    def test_default_tarball_path_uses_dataset_id(self):
        import geardb

        path = geardb.Dataset(id="XYZ", dtype="single-cell-rnaseq").get_tarball_path()
        assert path.endswith("/www/datasets/XYZ.tar.gz")
        spatial = geardb.Dataset(id="XYZ", dtype="spatial").get_tarball_path()
        assert spatial.endswith("/www/datasets/spatial/XYZ.tar.gz")
