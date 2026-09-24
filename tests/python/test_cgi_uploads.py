"""Dataset upload CGIs: file-type checks, complete-file checks and upload status reporting."""

import json
from pathlib import Path

import pytest

from helpers.cgi_harness import run_cgi

FAKES_DIR = Path(__file__).resolve().parent / "fakes"
SPATIAL_FAKE = {"gear.spatialhandler": FAKES_DIR / "uploads_spatialhandler.py"}
SHARE_UID = "pytestshare"


class TestStoreExpressionDataset:
    """store_expression_dataset.cgi saves the uploaded file as <share_uid>.<ext> in the upload area."""

    def run(self, session_id, filename, content=b"dataset bytes", dataset_format="mex_3tab",
            spatial_format=None, expected_size=None, logged_in=True):
        form = {"session_id": session_id, "share_uid": SHARE_UID, "dataset_format": dataset_format,
                "expected_size": len(content) if expected_size is None else expected_size}
        if spatial_format:
            form["spatial_format"] = spatial_format
        db = {"sessions": {session_id: 1}} if logged_in else {}
        return run_cgi("store_expression_dataset.cgi", form=form, files={"dataset_file": (filename, content)},
                       db=db, fake_modules=SPATIAL_FAKE)

    @pytest.fixture
    def share_dir(self, upload_session):
        session_id, path = upload_session
        (path / SHARE_UID).mkdir()
        return session_id, path / SHARE_UID

    @pytest.mark.parametrize("dataset_format, filename, saved_ext, spatial_format", [
        ("mex_3tab", "example.tar", "tar", None),
        ("mex_3tab", "bundle.TGZ", "tar.gz", None),
        ("mex_3tab", "bundle.tar.gz", "tar.gz", None),
        ("mex_3tab", "bundle.zip", "zip", None),
        ("excel", "data.xlsx", "xlsx", None),
        ("excel", "DATA.XLSX", "xlsx", None),
        ("rds", "seurat.RDS", "rds", None),
        ("h5ad", "adata.h5ad", "h5ad", None),
        ("spatial", "visium.tar", "tar", "visium"),
        ("spatial", "visium.tar.gz", "tar.gz", "visium"),
    ])
    def test_accepted_files_are_saved_under_a_fixed_name(self, share_dir, dataset_format, filename,
                                                         saved_ext, spatial_format):
        session_id, path = share_dir
        content = b"not really a " + filename.encode()
        result = self.run(session_id, filename, content, dataset_format=dataset_format,
                          spatial_format=spatial_format)
        assert result.json() == {"success": 1, "message": "Dataset file saved successfully."}
        saved = path / f"{SHARE_UID}.{saved_ext}"
        assert saved.read_bytes() == content
        assert [p.name for p in path.iterdir() if p.name != "status.json"] == [saved.name]
        status = json.loads((path / "status.json").read_text())
        assert status["status"] == "uploaded" and status["progress"] == 0

    def test_example_mex_tar_is_stored_intact(self, share_dir, example_mex_tar):
        session_id, path = share_dir
        content = example_mex_tar.read_bytes()
        result = self.run(session_id, "example_mex.tar", content)
        assert result.json()["success"] == 1
        assert (path / f"{SHARE_UID}.tar").read_bytes() == content

    @pytest.mark.parametrize("dataset_format, filename, message_fragment, spatial_format", [
        ("mex_3tab", "bundle.rar", ".tar, .tar.gz or .zip", None),
        ("excel", "old.xls", "Legacy .xls files are not supported", None),
        ("excel", "data.csv", "Expected .xlsx", None),
        ("h5ad", "x.rds", "Expected .h5ad", None),
        ("rds", "x.h5ad", "Expected .rds", None),
        ("spatial", "visium.zip", "Expected .tar or .tar.gz", "visium"),
        ("spatial", "visium.tar", "Invalid spatial format", "not-a-platform"),
    ])
    def test_wrong_file_types_are_rejected(self, share_dir, dataset_format, filename, message_fragment,
                                           spatial_format):
        session_id, path = share_dir
        result = self.run(session_id, filename, dataset_format=dataset_format, spatial_format=spatial_format)
        body = result.json()
        assert body["success"] == 0 and message_fragment in body["message"]
        assert list(path.iterdir()) == []

    def test_truncated_upload_is_rejected_and_removed(self, share_dir):
        session_id, path = share_dir
        result = self.run(session_id, "bundle.tar.gz", b"12345", expected_size=999)
        body = result.json()
        assert body["success"] == 0
        assert "incomplete" in body["message"] and "expected 999 bytes, received 5" in body["message"]
        assert list(path.iterdir()) == []

    def test_anonymous_upload_is_refused(self, share_dir):
        session_id, path = share_dir
        result = self.run(session_id, "bundle.tar.gz", logged_in=False)
        assert result.json() == {"success": 0, "message": "Only logged in users can upload datasets."}
        assert list(path.iterdir()) == []

    def test_missing_share_uid(self, upload_session):
        session_id, _ = upload_session
        result = run_cgi("store_expression_dataset.cgi", db={"sessions": {session_id: 1}},
                         form={"session_id": session_id, "dataset_format": "h5ad"},
                         files={"dataset_file": ("x.h5ad", b"data")})
        body = result.json()
        assert body["success"] == 0 and "share_uid missing" in body["message"]


class TestCheckDatasetProcessingStatus:
    """check_dataset_processing_status.cgi reads www/uploads/files/<session>/<share>/status.json."""

    def run(self, **query):
        return run_cgi("check_dataset_processing_status.cgi", query=query)

    def write_status(self, upload_session, status):
        session_id, path = upload_session
        (path / SHARE_UID).mkdir()
        (path / SHARE_UID / "status.json").write_text(json.dumps(status))
        return session_id

    def test_complete_reports_full_progress(self, upload_session):
        session_id = self.write_status(upload_session, {"job_id": "j", "status": "complete", "progress": 40})
        body = self.run(session_id=session_id, share_uid=SHARE_UID).json()
        assert body["status"] == "complete" and body["progress"] == 100

    def test_processing_job_with_uuid_is_returned_unchanged(self, upload_session):
        status = {"job_id": "0b8e8b3a-4f7c-4d4a-9a55-2f1c7f0b9d11", "status": "processing",
                  "message": "Processing the dataset.", "progress": 30}
        session_id = self.write_status(upload_session, status)
        result = self.run(session_id=session_id, share_uid=SHARE_UID)
        assert result.json() == status
        assert "TypeError" not in result.stderr

    def test_empty_status_is_not_treated_as_complete(self, upload_session):
        session_id = self.write_status(upload_session, {"job_id": None, "status": "", "progress": 0})
        body = self.run(session_id=session_id, share_uid=SHARE_UID).json()
        assert body["status"] == "" and body["progress"] == 0

    def test_processing_with_dead_process_is_an_error(self, upload_session):
        status = {"job_id": None, "process_id": 999999, "status": "processing", "progress": 10}
        session_id = self.write_status(upload_session, status)
        body = self.run(session_id=session_id, share_uid=SHARE_UID).json()
        assert body["status"] == "error" and body["progress"] == 0

    def test_missing_share_uid(self, upload_session):
        session_id, _ = upload_session
        body = self.run(session_id=session_id).json()
        assert body["status"] == "error" and body["message"] == "No share_uid provided."

    def test_missing_session_id(self):
        body = self.run(share_uid=SHARE_UID).json()
        assert body["status"] == "error" and body["message"] == "No session_id provided."

    def test_no_status_file(self, upload_session):
        session_id, _ = upload_session
        body = self.run(session_id=session_id, share_uid=SHARE_UID).json()
        assert body["status"] == "error" and "No status file found" in body["message"]

    def test_path_outside_uploads_is_refused(self):
        body = self.run(session_id="../..", share_uid="cgi").json()
        assert body["status"] == "error" and body["message"] == "Invalid session_id or share_uid."


class TestGetUploadsInProgress:
    def test_skips_directories_without_metadata(self, upload_session):
        session_id, path = upload_session
        good = path / "goodshare"
        good.mkdir()
        (good / "metadata.json").write_text(json.dumps(
            {"dataset_uid": "DS-1", "dataset_type": "single-cell-rnaseq", "title": "My upload"}))
        (path / "leftover").mkdir()     # e.g. partly recreated by a worker after the upload was deleted

        result = run_cgi("get_uploads_in_progress.cgi", query={"session_id": session_id})
        body = result.json()
        assert body["success"] == 1
        (upload,) = body["uploads"]
        assert upload["share_id"] == "goodshare"
        assert upload["dataset_id"] == "DS-1"
        assert upload["title"] == "My upload"
        assert upload["status"] == "metadata uploaded" and upload["load_step"] == "upload-dataset"

    def test_datafile_and_status_move_the_upload_along(self, upload_session):
        session_id, path = upload_session
        share = path / "share2"
        share.mkdir()
        (share / "metadata.json").write_text(json.dumps({"dataset_uid": "DS-2", "dataset_type": "bulk-rnaseq"}))
        (share / "share2.tar.gz").write_bytes(b"x")
        (share / "status.json").write_text(json.dumps({"status": "complete"}))
        (upload,) = run_cgi("get_uploads_in_progress.cgi", query={"session_id": session_id}).json()["uploads"]
        assert upload["status"] == "processed" and upload["load_step"] == "finalize-dataset"

    def test_no_upload_directory(self):
        body = run_cgi("get_uploads_in_progress.cgi", query={"session_id": "pytest-no-such-session"}).json()
        assert body == {"success": 1, "uploads": [], "message": ""}

    def test_missing_session_id(self):
        body = run_cgi("get_uploads_in_progress.cgi").json()
        assert body["success"] == 0 and body["message"] == "No session_id provided"
