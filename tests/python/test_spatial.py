"""
In-process tests for the spatial upload handlers (gear.spatialhandler) and
gear.spatial_processor.process_spatial_synchronously.

Needs the SpatialData stack (pip install -r requirements-spatial.txt); skipped otherwise.
"""

import io
import json
import shutil
import tarfile
from pathlib import Path

import pytest

pytest.importorskip("spatialdata")
pytest.importorskip("spatialdata_io")

import gear.spatial_processor as sp  # noqa: E402
from gear.spatialhandler import SPATIALTYPE2CLASS  # noqa: E402
from gear.utils.job_coordination import UploadCancelledError  # noqa: E402

pytestmark = pytest.mark.spatial

WRONG_TYPE_MESSAGE = "must be a .tar or .tar.gz file."


def tar_bytes(members: dict[str, bytes], gzipped: bool = False) -> bytes:
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz" if gzipped else "w") as tf:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    return buffer.getvalue()


README_TAR = tar_bytes({"readme.txt": b"not a real spatial dataset\n"})
README_TAR_GZ = tar_bytes({"readme.txt": b"not a real spatial dataset\n"}, gzipped=True)

# (filename, contents): every combination of real compression and extension
TARBALL_VARIANTS = {
    "plain_tar": ("upload.tar", README_TAR),
    "gzipped_tar_gz": ("upload.tar.gz", README_TAR_GZ),
    "plain_tar_named_tar_gz": ("upload.tar.gz", README_TAR),
    "gzipped_tar_named_tar": ("upload.tar", README_TAR_GZ),
}


@pytest.fixture
def extract_dir(tmp_path) -> Path:
    # Handlers extract into <extract_dir>/files and treat <extract_dir> as the upload's staging
    #  directory, stopping (UploadCancelledError) if it does not exist
    path = tmp_path / "ex"
    path.mkdir()
    return path


@pytest.mark.parametrize("variant", sorted(TARBALL_VARIANTS))
@pytest.mark.parametrize("spatial_type", sorted(SPATIALTYPE2CLASS))
def test_handler_opens_tarball_regardless_of_extension(spatial_type, variant, tmp_path, extract_dir):
    filename, contents = TARBALL_VARIANTS[variant]
    path = tmp_path / filename
    path.write_bytes(contents)

    handler = SPATIALTYPE2CLASS[spatial_type]()
    with pytest.raises(Exception) as excinfo:
        handler.process_file(str(path), extract_dir=str(extract_dir), organism_id=1)

    # The archive opened and was read; it failed later for lack of the platform's required files
    assert not isinstance(excinfo.value, tarfile.ReadError), excinfo.value
    assert not isinstance(excinfo.value, UploadCancelledError), excinfo.value
    assert WRONG_TYPE_MESSAGE not in str(excinfo.value)
    assert (extract_dir / "files" / "readme.txt").is_file()


@pytest.mark.parametrize("spatial_type", sorted(SPATIALTYPE2CLASS))
def test_handler_rejects_zip(spatial_type, tmp_path, extract_dir):
    path = tmp_path / "upload.zip"
    path.write_bytes(b"PK\x05\x06" + b"\x00" * 18)  # an empty zip archive

    handler = SPATIALTYPE2CLASS[spatial_type]()
    with pytest.raises(Exception, match=WRONG_TYPE_MESSAGE.replace(".", r"\.")):
        handler.process_file(str(path), extract_dir=str(extract_dir), organism_id=1)


@pytest.mark.parametrize("spatial_type, expected_names", [
    ("cosmx", ["exprMat_file.csv", "metadata_file.csv", "fov_positions_file.csv", "'CellComposite' directory", "'CellLabels' directory"]),
    ("curio", ["anndata.h5ad", "cluster_assignment.txt", "Metrics.csv", "variable_features_clusters.txt", "variable_features_spatial_moransi.txt"]),
    ("visium", ["filtered_feature_bc_matrix.h5", "clusters.csv", "spatial/tissue_hires_image.png", "spatial/scalefactors_json.json", "spatial/tissue_positions.csv"]),
    ("visium_hd", ["clusters.csv", "binned_outputs/feature_slice.h5", "binned_outputs/square_008um/filtered_feature_bc_matrix.h5", "binned_outputs/square_008um/spatial/tissue_hires_image.png"]),
    ("xenium", ["experiment.xenium", "cell_feature_matrix.h5", "cells.parquet", "morphology_focus.ome.tif", "clusters.csv"]),
])
def test_handler_reports_every_missing_input(spatial_type, expected_names, tmp_path, extract_dir):
    path = tmp_path / "upload.tar"
    path.write_bytes(README_TAR)

    handler = SPATIALTYPE2CLASS[spatial_type]()
    with pytest.raises(Exception, match="^Required .* files not found in tarball") as excinfo:
        handler.process_file(str(path), extract_dir=str(extract_dir), organism_id=1)
    for name in expected_names:
        assert name in str(excinfo.value)


def test_missing_input_report_skips_inputs_that_are_present(tmp_path, extract_dir):
    path = tmp_path / "upload.tar"
    path.write_bytes(tar_bytes({
        "experiment.xenium": b"{}",
        "clusters.csv": b"Barcode,Cluster\n",
        "morphology_focus/morphology_focus_0000.ome.tif": b"fake tif bytes",
    }))

    handler = SPATIALTYPE2CLASS["xenium"]()
    with pytest.raises(Exception) as excinfo:
        handler.process_file(str(path), extract_dir=str(extract_dir), organism_id=1)
    message = str(excinfo.value)
    assert "cell_feature_matrix.h5" in message and "cells.parquet" in message
    assert "experiment.xenium" not in message
    assert "clusters.csv" not in message
    assert "morphology_focus" not in message


@pytest.mark.parametrize("counts_name", ["filtered_feature_bc_matrix.h5", "Run1_filtered_feature_bc_matrix.h5"])
def test_visium_counts_file_prefix_is_optional(counts_name, tmp_path, extract_dir):
    path = tmp_path / "upload.tar"
    path.write_bytes(tar_bytes({counts_name: b"not really hdf5", "clusters.csv": b"Barcode,Cluster\n"}))

    handler = SPATIALTYPE2CLASS["visium"]()
    with pytest.raises(Exception, match="^Required Visium files not found") as excinfo:
        handler.process_file(str(path), extract_dir=str(extract_dir), organism_id=1)
    assert "filtered_feature_bc_matrix.h5" not in str(excinfo.value)
    assert "spatial/scalefactors_json.json" in str(excinfo.value)


def test_gzipped_readme_tar_is_really_gzipped():
    assert README_TAR_GZ[:2] == b"\x1f\x8b"
    assert README_TAR[:2] != b"\x1f\x8b"


# ---------------------------------------------------------------------------
# process_spatial_synchronously: cancellation when the upload is deleted
# ---------------------------------------------------------------------------

@pytest.fixture
def spatial_staging(tmp_path, monkeypatch) -> Path:
    monkeypatch.setattr(sp.geardb, "get_organism_id_by_taxon_id", lambda taxid: 1)
    path = tmp_path / "sess" / "SHARE"
    path.mkdir(parents=True)
    (path / "metadata.json").write_text(json.dumps({"sample_taxid": 10090}))
    (path / "SHARE.tar").write_bytes(
        tar_bytes(
            {
                "clusters.csv": b"Barcode,Cluster\nAAAC-1,1\n",
                "filtered_feature_bc_matrix.h5": b"not really hdf5",
                "spatial/tissue_positions.csv": b"barcode,in_tissue\n",
            }
        )
    )
    return path


def test_process_spatial_cancelled_when_upload_deleted_during_extraction(spatial_staging, monkeypatch):
    original_extract = tarfile.TarFile.extract
    calls = []

    def extract_then_delete(self, *args, **kwargs):
        value = original_extract(self, *args, **kwargs)
        calls.append(args)
        if len(calls) == 1:
            shutil.rmtree(spatial_staging)
        return value

    monkeypatch.setattr(tarfile.TarFile, "extract", extract_then_delete)

    result = sp.process_spatial_synchronously(
        job_id="job-1",
        share_uid="SHARE",
        staging_area=spatial_staging,
        status_file=spatial_staging / "status.json",
        spatial_format="visium",
        perform_primary_analysis=False,
    )

    assert result["success"] == 0, result
    assert result.get("cancelled") is True, result
    assert len(calls) == 1, "extraction continued after the upload was deleted"
    assert not spatial_staging.exists(), "the deleted staging directory was recreated"


def test_process_spatial_missing_metadata(tmp_path):
    staging = tmp_path / "sess" / "SHARE"
    staging.mkdir(parents=True)
    result = sp.process_spatial_synchronously(
        "job-1", "SHARE", staging, staging / "status.json", "visium", False
    )
    assert result["success"] == 0
    assert "No metadata JSON file found" in result["message"]


def test_cosmx_decompresses_gzipped_members(tmp_path, extract_dir):
    import gzip

    counts = b"fov,cell_ID,GeneA\n1,1,3\n"
    labels = b"fake tif bytes"
    path = tmp_path / "upload.tar.gz"
    path.write_bytes(tar_bytes({
        "RUN42_exprMat_file.csv.gz": gzip.compress(counts),
        "RUN42-CellLabels/CellLabels_F001.tif.gz": gzip.compress(labels),
        "notes.txt": b"plain member\n",
    }, gzipped=True))

    handler = SPATIALTYPE2CLASS["cosmx"]()
    with pytest.raises(Exception) as excinfo:
        handler.process_file(str(path), extract_dir=str(extract_dir), organism_id=1)
    # Fails later for lack of the other CosMx files, not while extracting
    assert not isinstance(excinfo.value, (tarfile.ReadError, UploadCancelledError)), excinfo.value

    files = extract_dir / "files"
    # The run's prefix is replaced by STANDARD_DATASET_ID ("spatialdata")
    assert (files / "spatialdata_exprMat_file.csv").read_bytes() == counts
    assert (files / "CellLabels" / "CellLabels_F001.tif").read_bytes() == labels
    assert (files / "notes.txt").is_file()
    assert not list(files.rglob("*.gz"))


COSMX_MEMBERS = {
    "exprMat_file.csv": b"fov,cell_ID,GeneA\n1,1,3\n",
    "metadata_file.csv": b"fov,cell_ID\n1,1\n",
    "fov_positions_file.csv": b"fov,x_global_px,y_global_px\n1,0,0\n",
    "CellComposite/CellComposite_F001.jpg": b"fake jpg bytes",
    "CellLabels/CellLabels_F001.tif": b"fake tif bytes",
}


@pytest.mark.parametrize("top_dir", ["", "Run1/"])
def test_cosmx_accepts_unprefixed_files(tmp_path, extract_dir, top_dir):
    path = tmp_path / "upload.tar.gz"
    path.write_bytes(tar_bytes({top_dir + name: data for name, data in COSMX_MEMBERS.items()}, gzipped=True))

    handler = SPATIALTYPE2CLASS["cosmx"]()
    with pytest.raises(Exception) as excinfo:
        handler.process_file(str(path), extract_dir=str(extract_dir), organism_id=1)
    # Fails later on the fake metadata, not because required files or directories were not found
    assert "RNA_Analysis_Neighborhood" in str(excinfo.value), excinfo.value

    files = extract_dir / "files"
    for suffix in ["exprMat_file.csv", "metadata_file.csv", "fov_positions_file.csv"]:
        assert (files / f"spatialdata_{suffix}").is_file()
    assert (files / "CellComposite" / "CellComposite_F001.jpg").is_file()
    assert (files / "CellLabels" / "CellLabels_F001.tif").is_file()


def test_cosmx_missing_file_error_uses_upload_names(tmp_path, extract_dir):
    members = {name: data for name, data in COSMX_MEMBERS.items() if name != "metadata_file.csv"}
    path = tmp_path / "upload.tar.gz"
    path.write_bytes(tar_bytes(members, gzipped=True))

    handler = SPATIALTYPE2CLASS["cosmx"]()
    with pytest.raises(Exception) as excinfo:
        handler.process_file(str(path), extract_dir=str(extract_dir), organism_id=1)
    message = str(excinfo.value)
    assert "metadata_file.csv" in message
    assert "spatialdata_" not in message



@pytest.mark.parametrize("extension, gzipped", [("tar.gz", True), ("tar", False)])
def test_truncated_archive_is_reported_as_incomplete(spatial_staging, extension, gzipped):
    """A spatial archive cut short (#1033) explains the likely cause instead of a raw Python error."""
    members = {
        "clusters.csv": b"Barcode,Cluster\n" + b"AAAC-1,1\n" * 20000,
        "filtered_feature_bc_matrix.h5": b"not really hdf5" * 5000,
    }
    (spatial_staging / "SHARE.tar").unlink()
    full = tar_bytes(members, gzipped=gzipped)
    (spatial_staging / f"SHARE.{extension}").write_bytes(full[: int(len(full) * 0.5)])

    result = sp.process_spatial_synchronously(
        job_id="job-1",
        share_uid="SHARE",
        staging_area=spatial_staging,
        status_file=spatial_staging / "status.json",
        spatial_format="visium",
        perform_primary_analysis=False,
    )
    assert result["success"] == 0, result
    assert not result.get("cancelled")
    assert result["message"].startswith("The uploaded archive appears to be incomplete"), result["message"]
    status = json.loads((spatial_staging / "status.json").read_text())
    assert status["status"] == "error" and status["message"] == result["message"]
