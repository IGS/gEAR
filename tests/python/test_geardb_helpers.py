"""lib/geardb.py helpers that don't need a database."""

import importlib

import pytest


def test_geardb_imports_without_mysql():
    """conftest.py replaces MySQLDB.connect, so importing geardb needs neither MySQL nor gear.ini."""
    geardb = importlib.import_module("geardb")
    assert hasattr(geardb, "Dataset")
    assert callable(geardb.get_dataset_by_id)


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
