"""gear.analysis.ZarrAdapter path guard, and exception messages that never contain filesystem paths."""

import os
from pathlib import Path

import pytest

import gear.analysis as analysis


@pytest.fixture
def www(tmp_path, monkeypatch) -> Path:
    """
    A www/ whose datasets directory is a symlink to a separate "mounted disk", as on the hosts
    with hyperdisks, plus a normal (non-symlinked) uploads directory.
    """
    disk = tmp_path / "mnt" / "datastore" / "datasets"
    (disk / "spatial" / "DS1.zarr" / "tables").mkdir(parents=True)
    www = tmp_path / "gEAR" / "www"
    (www / "uploads" / "files" / "sess" / "SHARE" / "SHARE.zarr").mkdir(parents=True)
    os.symlink(disk, www / "datasets")
    monkeypatch.setattr(analysis, "ZARR_ALLOWED_DIRS", (
        www / "datasets", www / "datasets" / "spatial", www / "uploads" / "files"))
    return www


def test_dataset_on_a_symlinked_disk_is_allowed(www):
    path = www / "datasets" / "spatial" / "DS1.zarr"
    assert path.resolve().is_relative_to(www.parent.parent / "mnt")   # really lives outside www
    assert analysis.ZarrAdapter(path).zarr_path == path


def test_upload_staging_area_is_allowed(www):
    path = www / "uploads" / "files" / "sess" / "SHARE" / "SHARE.zarr"
    assert analysis.ZarrAdapter(path).zarr_path == path


@pytest.mark.parametrize("relative", ["datasets/spatial/../../../../etc", "analyses/DS1.zarr", "."])
def test_paths_outside_the_allowed_directories_are_refused(www, relative, capsys):
    with pytest.raises(ValueError) as excinfo:
        analysis.ZarrAdapter(www / relative)
    assert str(excinfo.value) == "The requested dataset is not in an allowed location."
    assert "/" not in str(excinfo.value)
    # The path is still logged for whoever debugs it
    assert "refusing zarr path" in capsys.readouterr().err


def test_symlink_inside_an_allowed_directory_cannot_escape(www, tmp_path):
    secret = tmp_path / "secret.zarr"
    secret.mkdir()
    os.symlink(secret, www / "datasets" / "spatial" / "evil.zarr")
    with pytest.raises(ValueError):
        analysis.ZarrAdapter(www / "datasets" / "spatial" / "evil.zarr")


# get_sdata() imports spatialdata, so only this test needs the spatial stack. The mark lets
#  "pytest -m spatial" select it; importorskip skips it (rather than failing) where the stack
#  isn't installed, e.g. the core CI job. The rest of this module runs everywhere.
@pytest.mark.spatial
def test_missing_store_message_has_no_path(www):
    pytest.importorskip("spatialdata")
    missing = analysis.ZarrAdapter(www / "datasets" / "spatial" / "NOPE.zarr")
    with pytest.raises(FileNotFoundError) as excinfo:
        missing.get_sdata()
    assert str(excinfo.value) == "Dataset NOPE was not found"


def test_missing_table_message_has_no_path(www):
    no_table = analysis.ZarrAdapter(www / "datasets" / "spatial" / "DS1.zarr")
    with pytest.raises(FileNotFoundError) as excinfo:
        no_table.get_adata()
    assert str(www) not in str(excinfo.value) and "/" not in str(excinfo.value)


def test_missing_primary_dataset_message_has_no_path(monkeypatch, tmp_path):
    class FakeDataset:
        def __init__(self, **kwargs):
            pass

        def get_file_path(self):
            return str(tmp_path / "datasets" / "DS9.h5ad")

    monkeypatch.setattr(analysis, "Dataset", FakeDataset)
    with pytest.raises(FileNotFoundError) as excinfo:
        analysis.get_primary_analysis("DS9")
    assert str(excinfo.value) == "No h5 file found for dataset DS9"
