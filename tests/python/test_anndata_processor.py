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
import sys
import tarfile
import zipfile
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
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
    assert isinstance(metadata["primary_variable_candidates"], dict)


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


def test_threetab_ontology_names_added(staging, example_3tab_tar_gz, monkeypatch):
    """Ontology term IDs in the observations get a name column, and are listed in metadata.json (#422)."""
    from gear import ontology
    monkeypatch.setattr(ontology, "_read_cache", lambda term_ids: {})
    monkeypatch.setattr(ontology, "_write_cache", lambda labels: None)
    monkeypatch.setattr(ontology, "_fetch_from_ols", lambda term_ids: {"UBERON:0002227": "spiral organ of cochlea"})
    monkeypatch.setattr(ontology, "get_recommended_obs_columns", lambda: ["tissue_ontology", "age_value"])

    with tarfile.open(example_3tab_tar_gz) as tf:
        members = {m.name: tf.extractfile(m).read() for m in tf.getmembers() if m.isfile()}
    lines = members["observations.tab"].decode().splitlines()
    lines = [lines[0] + "\ttissue_ontology"] + [line + "\tUBERON:0002227" for line in lines[1:]]
    members["observations.tab"] = ("\n".join(lines) + "\n").encode()
    write_tar(staging / "SHARE.tar.gz", members, gzipped=True)

    assert_success(run(staging), staging)
    adata = read_h5ad(staging)
    assert set(adata.obs["tissue_ontology_name"]) == {"spiral organ of cochlea"}
    metadata = json.loads((staging / "metadata.json").read_text())
    assert metadata["ontology_terms"] == {"tissue_type": ["UBERON:0002227"]}
    assert len(metadata["sample_metadata_warnings"]) == 1 and "age_value" in metadata["sample_metadata_warnings"][0]


# ---------------------------------------------------------------------------
# Bad or missing archives
# ---------------------------------------------------------------------------

def test_not_an_archive(staging):
    (staging / "SHARE.tar").write_bytes(b"hello")
    result = run(staging)
    assert result["success"] == 0
    assert not result.get("cancelled")
    assert result["message"].startswith("The uploaded archive could not be read")
    assert "tar -tf" in result["message"]
    status = json.loads((staging / "status.json").read_text())
    assert status["status"] == "error"



def truncated(data: bytes, fraction: float) -> bytes:
    return data[: int(len(data) * fraction)]


@pytest.mark.parametrize("fraction", [0.3, 0.9, 0.99])
@pytest.mark.parametrize("extension, archive", [
    ("tar.gz", lambda members, path: write_tar(path, members, gzipped=True)),
    ("tar", lambda members, path: write_tar(path, members)),
    ("zip", lambda members, path: write_zip(path, members)),
])
def test_truncated_archive_is_reported_as_incomplete(staging, tmp_path, mex_members, extension, archive, fraction):
    """An archive cut short (#1032, #1033) explains the likely cause instead of an internal error."""
    full = archive(mex_members, tmp_path / f"full.{extension}").read_bytes()
    (staging / f"SHARE.{extension}").write_bytes(truncated(full, fraction))

    result = run(staging)
    assert result["success"] == 0
    assert not result.get("cancelled")
    message = result["message"]
    assert "internal error" not in message
    assert message.startswith("The uploaded archive")
    assert ("unzip -t" if extension == "zip" else "tar -tf") in message
    status = json.loads((staging / "status.json").read_text())
    assert status["status"] == "error" and status["message"] == message


def test_truncated_tar_gz_says_incomplete(staging, tmp_path, mex_members):
    full = write_tar(tmp_path / "full.tar.gz", mex_members, gzipped=True).read_bytes()
    (staging / "SHARE.tar.gz").write_bytes(truncated(full, 0.9))
    assert run(staging)["message"].startswith("The uploaded archive appears to be incomplete")

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


# ---------------------------------------------------------------------------
# RDS / Seurat uploads
# ---------------------------------------------------------------------------
# rpy2 and R aren't installed for the tests, so gear.seuratuploader is imported with rpy2
# stubbed out. Only the R step (seurat_to_anndata) and the network gene mapping are replaced;
# the processor and layer_to_X run as in production.

@pytest.fixture
def seurat_uploader(monkeypatch):
    """gear.seuratuploader imported with its rpy2 imports stubbed."""
    from unittest import mock

    for name in (
        "rpy2",
        "rpy2.rinterface_lib",
        "rpy2.rinterface_lib.callbacks",
        "rpy2.robjects",
        "rpy2.robjects.packages",
    ):
        monkeypatch.setitem(sys.modules, name, mock.MagicMock())
    monkeypatch.delitem(sys.modules, "gear.seuratuploader", raising=False)
    import gear.seuratuploader as seurat_uploader

    monkeypatch.setattr(seurat_uploader, "genes_to_ensembl", lambda adata, taxid=None: adata)
    return seurat_uploader


def converted_seurat(layers: dict[str, np.ndarray]) -> anndata.AnnData:
    """What anndataR writes for a Seurat object: no X, and each assay layer under its own name."""
    obs = pd.DataFrame({"cluster": ["a", "b", "a"]}, index=["c1", "c2", "c3"])
    var = pd.DataFrame(index=["Xkr4", "Sox17"])
    return anndata.AnnData(obs=obs, var=var, layers=layers)


def run_seurat(staging: Path, seurat_uploader, monkeypatch, layers: dict[str, np.ndarray]) -> dict:
    (staging / f"{SHARE_UID}.rds").write_bytes(b"not read: the R step is replaced")

    def fake_seurat_to_anndata(file_path, share_name, output_dir=".", progress_callback=None):
        output_path = Path(output_dir) / f"tmp_{share_name}.h5ad"
        converted_seurat(layers).write_h5ad(output_path)
        return str(output_path)

    monkeypatch.setattr(seurat_uploader, "seurat_to_anndata", fake_seurat_to_anndata)
    return make_processor(staging).process("rds")


COUNTS = np.array([[1, 0], [3, 2], [0, 5]], dtype=np.float32)
LOG_NORMALIZED = np.log1p(COUNTS)


def test_seurat_data_layer_becomes_x(staging, seurat_uploader, monkeypatch):
    """The normalized "data" layer is used for expression values, never the raw counts."""
    result = run_seurat(
        staging, seurat_uploader, monkeypatch, {"counts": COUNTS, "data": LOG_NORMALIZED}
    )

    assert result["success"] == 1, result
    adata = read_h5ad(staging)
    np.testing.assert_allclose(np.asarray(adata.X), LOG_NORMALIZED)


@pytest.mark.parametrize(
    "layers, advice",
    [
        # NormalizeData() was never run on the object
        ({"counts": COUNTS}, "NormalizeData"),
        # Seurat v5 objects split by sample (split() or integration) without JoinLayers()
        ({"counts.1": COUNTS, "data.1": LOG_NORMALIZED, "data.2": LOG_NORMALIZED}, "JoinLayers"),
    ],
    ids=["not-normalized", "split-layers"],
)
def test_seurat_without_data_layer_is_explained(staging, seurat_uploader, monkeypatch, layers, advice):
    """
    A Seurat object with no layer named "data" fails with a message that names the layers found
    and says how to fix the object, instead of "An unexpected internal error ... 'data'".
    """
    result = run_seurat(staging, seurat_uploader, monkeypatch, layers)

    assert result["success"] == 0, result
    message = result["message"]
    assert "unexpected internal error" not in message
    assert advice in message
    for layer_name in layers:
        assert layer_name in message

    status = json.loads((staging / "status.json").read_text())
    assert status["status"] == "error"
    assert status["message"] == message
    # The raw counts are never used as expression values in place of the missing "data" layer
    assert not (staging / f"{SHARE_UID}.h5ad").exists()
