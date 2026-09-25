"""
In-process tests for gear.anndata_processor.AnndataProcessor on MEX and 3-tab archive uploads.

Each test uses tmp_path as the uploads base directory, with a staging area at
tmp_path/sess/SHARE holding the uploaded archive (saved as SHARE.<extension>, as
store_expression_dataset.cgi does) and a metadata.json.
"""

import gzip
import io
import json
import shutil
import tarfile
import zipfile
from pathlib import Path

import anndata
import pytest

import gear.anndata_processor as ap

SHARE_UID = "SHARE"
MEX_FILES = ("matrix.mtx", "barcodes.tsv", "genes.tsv")


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------

@pytest.fixture
def staging(tmp_path, monkeypatch) -> Path:
    """A staging area (tmp_path/sess/SHARE) under a temporary uploads base directory."""
    monkeypatch.setattr(ap, "UPLOADS_BASE_DIR", tmp_path)
    path = tmp_path / "sess" / SHARE_UID
    path.mkdir(parents=True)
    (path / "metadata.json").write_text("{}")
    return path


@pytest.fixture(scope="module")
def mex_members() -> dict[str, bytes]:
    """The contents of the MEX example's matrix.mtx, barcodes.tsv and genes.tsv, keyed by name."""
    from conftest import USER_TEMPLATES

    with tarfile.open(USER_TEMPLATES / "example_mex.tar") as tf:
        return {name: tf.extractfile(name).read() for name in MEX_FILES}


def make_processor(staging: Path) -> ap.AnndataProcessor:
    return ap.AnndataProcessor("job-1", SHARE_UID, staging, staging / "status.json", "DS")


def run(staging: Path) -> dict:
    return make_processor(staging).process("mex_3tab")


def write_tar(path: Path, members: dict[str, bytes], gzipped: bool = False) -> Path:
    with tarfile.open(path, "w:gz" if gzipped else "w") as tf:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    return path


def write_zip(path: Path, members: dict[str, bytes]) -> Path:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, data in members.items():
            zf.writestr(name, data)
    return path


def read_h5ad(staging: Path) -> anndata.AnnData:
    return anndata.read_h5ad(staging / f"{SHARE_UID}.h5ad")


def assert_success(result: dict, staging: Path) -> None:
    assert result["success"] == 1, result
    assert not result.get("cancelled")
    status = json.loads((staging / "status.json").read_text())
    assert status["status"] == "complete"
    assert status["progress"] == 100
    metadata = json.loads((staging / "metadata.json").read_text())
    assert "questionable_obs_columns" in metadata
    assert "obs_dtype_reviewed" in metadata


def assert_mex_example(adata: anndata.AnnData) -> None:
    assert adata.shape == (150, 27998)
    assert adata.var.index[0] == "ENSMUSG00000051951"
    assert "gene_symbol" in adata.var.columns
    assert adata.var["gene_symbol"].iloc[0] == "Xkr4"


def assert_cancelled(result: dict, staging: Path) -> None:
    assert result["success"] == 0, result
    assert result.get("cancelled") is True, result
    assert "deleted" in result["message"]
    assert not staging.exists(), "the deleted staging directory was recreated"


# ---------------------------------------------------------------------------
# MEX uploads
# ---------------------------------------------------------------------------

def test_mex_plain_tar(staging, example_mex_tar):
    shutil.copy(example_mex_tar, staging / "SHARE.tar")
    result = run(staging)
    assert_success(result, staging)
    assert_mex_example(read_h5ad(staging))


def test_mex_files_nested_in_folder(staging, mex_members):
    write_tar(
        staging / "SHARE.tar",
        {f"filtered_feature_bc_matrix/{name}": data for name, data in mex_members.items()},
    )
    result = run(staging)
    assert_success(result, staging)
    assert_mex_example(read_h5ad(staging))


def test_mex_cellranger_v3_gzipped_files(staging, mex_members):
    features = "".join(
        f"{line}\tGene Expression\n"
        for line in mex_members["genes.tsv"].decode().splitlines()
        if line
    ).encode()
    write_tar(
        staging / "SHARE.tar.gz",
        {
            "filtered_feature_bc_matrix/matrix.mtx.gz": gzip.compress(mex_members["matrix.mtx"]),
            "filtered_feature_bc_matrix/barcodes.tsv.gz": gzip.compress(mex_members["barcodes.tsv"]),
            "filtered_feature_bc_matrix/features.tsv.gz": gzip.compress(features),
        },
        gzipped=True,
    )
    result = run(staging)
    assert_success(result, staging)
    adata = read_h5ad(staging)
    assert_mex_example(adata)
    assert "feature_types" in adata.var.columns
    assert (adata.var["feature_types"] == "Gene Expression").all()


def test_mex_zip(staging, mex_members):
    write_zip(staging / "SHARE.zip", mex_members)
    result = run(staging)
    assert_success(result, staging)
    assert_mex_example(read_h5ad(staging))


def test_mex_gzipped_tar_named_tar(staging, mex_members):
    write_tar(staging / "SHARE.tar", mex_members, gzipped=True)
    with open(staging / "SHARE.tar", "rb") as fh:
        assert fh.read(2) == b"\x1f\x8b"  # really gzipped
    result = run(staging)
    assert_success(result, staging)
    assert_mex_example(read_h5ad(staging))


def test_mex_plain_tar_named_tar_gz(staging, example_mex_tar):
    shutil.copy(example_mex_tar, staging / "SHARE.tar.gz")
    result = run(staging)
    assert_success(result, staging)
    assert_mex_example(read_h5ad(staging))


# ---------------------------------------------------------------------------
# 3-tab uploads
# ---------------------------------------------------------------------------

def test_threetab_tar_gz(staging, example_3tab_tar_gz):
    shutil.copy(example_3tab_tar_gz, staging / "SHARE.tar.gz")
    result = run(staging)
    assert_success(result, staging)
    assert read_h5ad(staging).shape == (24, 31407)


def test_threetab_decompressed_plain_tar(staging, example_3tab_tar_gz):
    (staging / "SHARE.tar").write_bytes(gzip.decompress(example_3tab_tar_gz.read_bytes()))
    assert tarfile.is_tarfile(staging / "SHARE.tar")
    result = run(staging)
    assert_success(result, staging)
    assert read_h5ad(staging).shape == (24, 31407)


# ---------------------------------------------------------------------------
# Bad or missing archives
# ---------------------------------------------------------------------------

def test_not_an_archive(staging):
    (staging / "SHARE.tar").write_bytes(b"hello")
    result = run(staging)
    assert result["success"] == 0
    assert not result.get("cancelled")
    assert "tar archive could not be read" in result["message"]
    status = json.loads((staging / "status.json").read_text())
    assert status["status"] == "error"


def test_missing_archive(staging):
    result = run(staging)
    assert result["success"] == 0
    assert not result.get("cancelled")
    assert "could not be found" in result["message"]


def test_staging_area_outside_uploads_dir_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(ap, "UPLOADS_BASE_DIR", tmp_path / "uploads")
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    with pytest.raises(ap.ProcessingError, match="Invalid staging area path"):
        make_processor(outside)


# ---------------------------------------------------------------------------
# Cancellation: the user deletes the upload (staging directory) mid-processing
# ---------------------------------------------------------------------------

def _delete_after_first_call(monkeypatch, owner, attr, staging: Path) -> list:
    """Wrap owner.attr so it runs normally, then removes the staging directory after the first call."""
    original = getattr(owner, attr)
    calls = []

    def wrapper(*args, **kwargs):
        value = original(*args, **kwargs)
        calls.append(args)
        if len(calls) == 1:
            shutil.rmtree(staging)
        return value

    monkeypatch.setattr(owner, attr, wrapper)
    return calls


def test_cancel_during_tar_extraction(staging, example_mex_tar, monkeypatch):
    shutil.copy(example_mex_tar, staging / "SHARE.tar")
    calls = _delete_after_first_call(monkeypatch, tarfile.TarFile, "extract", staging)
    result = run(staging)
    assert_cancelled(result, staging)
    assert len(calls) == 1, "extraction continued after the upload was deleted"


def test_cancel_during_zip_extraction(staging, mex_members, monkeypatch):
    write_zip(staging / "SHARE.zip", mex_members)
    calls = _delete_after_first_call(monkeypatch, zipfile.ZipFile, "extract", staging)
    result = run(staging)
    assert_cancelled(result, staging)
    assert len(calls) == 1, "extraction continued after the upload was deleted"


def test_cancel_during_obs_sanitizing(staging, example_mex_tar, monkeypatch):
    shutil.copy(example_mex_tar, staging / "SHARE.tar")
    original = ap.standardize_and_sanitize_obs

    def delete_then_sanitize(obs):
        shutil.rmtree(staging)
        return original(obs)

    monkeypatch.setattr(ap, "standardize_and_sanitize_obs", delete_then_sanitize)
    result = run(staging)
    assert_cancelled(result, staging)


def test_status_updates_never_recreate_deleted_staging_dir(staging):
    processor = make_processor(staging)
    processor._update_status("processing", "started")
    assert (staging / "status.json").is_file()

    shutil.rmtree(staging)

    processor._update_status("error", "should not be written")
    processor._write_status_file()
    assert not staging.exists()

    with pytest.raises(ap.UploadCancelledError):
        processor._update_progress(50, "should not be written")
    assert not staging.exists()

    assert processor._cancelled_result()["cancelled"] is True
