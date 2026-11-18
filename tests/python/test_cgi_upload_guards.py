"""Upload CGIs that build paths from share_uid/session_id/dataset_uid: only plain names are accepted."""

import json
import shutil
import uuid
from pathlib import Path

import pytest

from helpers.cgi_harness import run_cgi

FAKES_DIR = Path(__file__).resolve().parent / "fakes"
METADATA_FAKE = {"gear.metadata": FAKES_DIR / "guards_metadata.py"}
WWW_DIR = Path(__file__).resolve().parents[2] / "www"
UPLOADS_DIR = WWW_DIR / "uploads" / "files"
DATASETS_DIR = WWW_DIR / "datasets"
SHARE_UID = "pytestshare"

BAD_IDS = ["..", ".", "../escape", "a/b", "x y", "/tmp/abs"]


@pytest.fixture
def share_dir(upload_session):
    session_id, path = upload_session
    share = path / SHARE_UID
    share.mkdir()
    (share / "metadata.json").write_text("{}")
    return session_id, share


class TestDeleteUploadInProgress:
    def run(self, session_id, share_uid, logged_in=True):
        db = {"sessions": {session_id: 1}} if logged_in else {}
        query = {"session_id": session_id}
        if share_uid is not None:
            query["share_uid"] = share_uid
        return run_cgi("delete_upload_in_progress.cgi", query=query, db=db)

    def test_owner_deletes_upload(self, share_dir):
        session_id, share = share_dir
        result = self.run(session_id, SHARE_UID)
        assert result.json() == {"success": 1, "message": "File deleted successfully."}
        assert not share.exists()
        assert share.parent.is_dir()

    @pytest.mark.parametrize("share_uid", BAD_IDS + [None, ""])
    def test_odd_share_uid_is_refused(self, share_dir, share_uid):
        session_id, share = share_dir
        result = self.run(session_id, share_uid)
        body = result.json()
        assert body["success"] == 0 and body["message"].startswith("Invalid share_uid")
        assert share.is_dir() and (share / "metadata.json").is_file()
        assert result.writes() == []

    def test_traversal_in_session_id_is_refused(self, share_dir):
        session_id, share = share_dir
        odd_session = f"../{session_id}"
        result = run_cgi("delete_upload_in_progress.cgi", db={"sessions": {odd_session: 1}},
                         query={"session_id": odd_session, "share_uid": SHARE_UID})
        assert result.json() == {"success": 0, "message": "Invalid session_id."}
        assert share.is_dir()

    def test_anonymous_is_refused(self, share_dir):
        session_id, share = share_dir
        result = self.run(session_id, SHARE_UID, logged_in=False)
        body = result.json()
        assert body["success"] == 0 and "log in" in body["message"]
        assert share.is_dir()
        assert result.writes() == []

    def test_missing_upload_directory(self, upload_session):
        session_id, _ = upload_session
        body = self.run(session_id, "nosuchshare").json()
        assert body["success"] == 0 and body["message"].startswith("Upload directory not found")


class TestStoreExpressionMetadata:
    def run(self, session_id, share_uid, logged_in=True, **fields):
        db = {"sessions": {session_id: 1}} if logged_in else {}
        form = {"session_id": session_id, "share_uid": share_uid, **fields}
        return run_cgi("store_expression_metadata.cgi", form=form, db=db)

    def test_metadata_is_saved(self, upload_session):
        session_id, path = upload_session
        result = self.run(session_id, SHARE_UID, dataset_uid="DS-1", title="My dataset", taxon_id="10090")
        assert result.json() == {"success": 1, "error": None}
        saved = json.loads((path / SHARE_UID / "metadata.json").read_text())
        assert saved["owner_id"] == 1
        assert saved["dataset_uid"] == "DS-1" and saved["share_uid"] == SHARE_UID
        assert saved["title"] == "My dataset" and saved["sample_taxid"] == "10090"
        assert result.writes() == []

    @pytest.mark.parametrize("share_uid", BAD_IDS)
    def test_odd_share_uid_is_refused(self, upload_session, share_uid):
        session_id, path = upload_session
        result = self.run(session_id, share_uid, dataset_uid="DS-1")
        assert result.json() == {"success": 0, "error": "Invalid share_uid"}
        assert list(path.iterdir()) == []
        assert result.writes() == []

    def test_odd_session_id_is_refused(self, upload_session):
        session_id, path = upload_session
        odd_session = f"{session_id}/../{session_id}"
        result = run_cgi("store_expression_metadata.cgi", db={"sessions": {odd_session: 1}},
                         form={"session_id": odd_session, "share_uid": SHARE_UID})
        assert result.json() == {"success": 0, "error": "Invalid session_id"}
        assert list(path.iterdir()) == []

    def test_anonymous_is_refused_without_creating_directories(self):
        session_id = f"pytest-{uuid.uuid4().hex[:12]}"
        try:
            result = self.run(session_id, SHARE_UID, logged_in=False)
            assert result.json() == {"success": 0, "error": "Invalid session_id"}
            assert not (UPLOADS_DIR / session_id).exists()
        finally:
            shutil.rmtree(UPLOADS_DIR / session_id, ignore_errors=True)

    def test_missing_session_id(self):
        result = run_cgi("store_expression_metadata.cgi", form={"share_uid": SHARE_UID})
        assert result.json() == {"success": 0, "error": "Invalid session_id"}


class TestFinalizeUploadedExpressionDataset:
    def run(self, session_id, share_uid=SHARE_UID, dataset_uid="DS-1", dataset_format="h5ad", **extra):
        form = {"session_id": session_id, "dataset_uid": dataset_uid, "dataset_format": dataset_format,
                "dataset_visibility": "private", **extra}
        if share_uid is not None:
            form["share_uid"] = share_uid
        return run_cgi("finalize_uploaded_expression_dataset.cgi", form=form,
                       db={"sessions": {session_id: 1}}, fake_modules=METADATA_FAKE)

    @pytest.fixture
    def staged(self, share_dir):
        session_id, share = share_dir
        (share / "metadata.json").write_text(json.dumps({"dataset_uid": "DS-1"}))
        (share / f"{SHARE_UID}.h5ad").write_bytes(b"h5ad bytes")
        (share / f"{SHARE_UID}.tar.gz").write_bytes(b"tarball bytes")
        (share / "status.json").write_text(json.dumps({"status": "complete"}))
        return session_id, share

    def assert_untouched(self, result, share):
        assert sorted(p.name for p in share.iterdir()) == sorted(
            ["metadata.json", f"{SHARE_UID}.h5ad", f"{SHARE_UID}.tar.gz", "status.json"])
        assert not result.logged("metadata.save_to_mysql")
        assert result.writes() == []

    def test_missing_share_uid(self, staged):
        session_id, share = staged
        result = self.run(session_id, share_uid=None)
        body = result.json()
        assert body["success"] == 0 and body["message"] == "Invalid share_uid: None"
        assert "Traceback" not in result.stderr
        self.assert_untouched(result, share)

    @pytest.mark.parametrize("share_uid", BAD_IDS + [f"../{SHARE_UID}", f"{SHARE_UID}/."])
    def test_odd_share_uid_is_refused(self, staged, share_uid):
        session_id, share = staged
        result = self.run(session_id, share_uid=share_uid)
        body = result.json()
        assert body["success"] == 0 and body["message"].startswith("Invalid share_uid")
        self.assert_untouched(result, share)

    @pytest.mark.parametrize("dataset_uid", ["..", "../escape", "a/b", "x y"])
    def test_odd_dataset_uid_is_refused(self, staged, dataset_uid):
        session_id, share = staged
        result = self.run(session_id, dataset_uid=dataset_uid)
        body = result.json()
        assert body["success"] == 0 and body["message"].startswith("Invalid dataset_uid")
        self.assert_untouched(result, share)

    def test_dataset_uid_cannot_move_files_out_of_datasets_dir(self, staged):
        # Without the guard the H5AD would land in www/uploads/files/<session>/escaped.h5ad
        session_id, share = staged
        result = self.run(session_id, dataset_uid=f"../uploads/files/{session_id}/escaped")
        assert result.json()["success"] == 0
        assert not (share.parent / "escaped.h5ad").exists()
        self.assert_untouched(result, share)

    def test_anonymous_is_refused(self, staged):
        session_id, share = staged
        result = run_cgi("finalize_uploaded_expression_dataset.cgi", fake_modules=METADATA_FAKE,
                         form={"session_id": session_id, "share_uid": SHARE_UID, "dataset_uid": "DS-1",
                               "dataset_format": "h5ad", "dataset_visibility": "private"})
        assert "log in" in result.json()["message"]
        self.assert_untouched(result, share)

    def test_staged_files_are_moved_into_the_datasets_dir(self, staged):
        session_id, share = staged
        (share / "analysis_pipeline.json").write_text("{}")
        dataset_uid = f"pytest-{uuid.uuid4().hex}"
        moved = [DATASETS_DIR / f"{dataset_uid}.{ext}" for ext in ("h5ad", "tar.gz", "pipeline.json")]
        try:
            result = self.run(session_id, dataset_uid=dataset_uid, dataset_format="mex_3tab",
                              perform_analysis_migration="1")
            assert result.json() == {"success": 1, "metadata_loaded": 1, "h5ad_migrated": 1,
                                     "userdata_migrated": 1, "primary_analysis_migrated": 1,
                                     "message": "All steps completed successfully."}
            assert moved[0].read_bytes() == b"h5ad bytes"
            assert moved[1].read_bytes() == b"tarball bytes"
            assert moved[2].read_text() == "{}"
            assert not share.exists()
            (saved,) = result.logged("metadata.save_to_mysql")
            assert saved["status"] == "completed" and saved["is_public"] == 0
        finally:
            for path in moved:
                path.unlink(missing_ok=True)
