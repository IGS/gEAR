"""POST /import/dataset/<share_uid>/status: auth, path safety and reading the upload's status.json."""

import json
import subprocess
import types

import pytest
from flask import Flask
from flask_restful import Api

import resources.dataset_processing as dataset_processing

SESSION = "SESS1"
SHARE = "a1b2c3d4"     # guid('short') in www/js/common.v2.js: the first block of a UUID


@pytest.fixture
def env(monkeypatch, tmp_path):
    """A test client with one valid session, uploads under tmp_path, and a recorded fake `ps`."""
    mod = dataset_processing
    ps_calls = []
    ps = types.SimpleNamespace(returncode=0)

    def fake_run(cmd, **kwargs):
        ps_calls.append(cmd)
        return types.SimpleNamespace(returncode=ps.returncode)

    uploads = tmp_path / "uploads"
    uploads.mkdir()
    monkeypatch.setattr(mod, "user_upload_file_base", uploads)
    monkeypatch.setattr(mod.geardb, "get_user_from_session_id",
                        lambda session_id=None: types.SimpleNamespace(id=1) if session_id == SESSION else None)
    monkeypatch.setattr(subprocess, "run", fake_run)

    app = Flask(__name__)
    Api(app).add_resource(mod.DatasetProcessingStatus, "/import/dataset/<share_uid>/status")
    client = app.test_client()
    client.set_cookie("gear_session_id", SESSION)
    return types.SimpleNamespace(client=client, tmp=tmp_path, uploads=uploads, ps=ps, ps_calls=ps_calls,
                                 module=mod, staging=uploads / SESSION / SHARE)


def post(env, share=SHARE, dataset_format="expression"):
    resp = env.client.post(f"/import/dataset/{share}/status", json={"dataset_format": dataset_format})
    return resp.status_code, resp.get_json()


def write_status(env, status, path=None):
    path = path or env.staging / "status.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(status))


ERROR_KEYS = {"success", "message", "status", "progress"}


def test_error_response_shape():
    body, status = dataset_processing._error_response("Nope", 418)
    assert status == 418
    assert body == {"success": False, "message": "Nope", "status": "error", "progress": 0}


def test_unknown_session(env):
    write_status(env, {"status": "complete"})
    env.client.set_cookie("gear_session_id", "NOT-A-SESSION")
    status, body = post(env)
    assert status == 401
    assert body == {"success": False, "message": "Invalid session", "status": "error", "progress": 0}


def test_missing_session_cookie(env):
    env.client.delete_cookie("gear_session_id")
    status, body = post(env)
    assert status == 401
    assert body["message"] == "Invalid session"


def test_session_that_is_not_a_plain_name(env, monkeypatch):
    # Even if the database accepted it, the session ID must not be usable to leave the uploads area
    monkeypatch.setattr(env.module.geardb, "get_user_from_session_id",
                        lambda session_id=None: types.SimpleNamespace(id=1))
    # "<uploads>/../x/status.json" would otherwise be read
    write_status(env, {"status": "complete", "secret": True}, env.tmp / "x" / "status.json")
    env.client.set_cookie("gear_session_id", "..")
    status, body = post(env, share="x")
    assert status == 401
    assert set(body) == ERROR_KEYS


@pytest.mark.parametrize("share", ["a.b", "%2E%2E", ".hidden", "bad share", "a;b"])
def test_share_id_that_is_not_a_plain_name(env, share):
    status, body = post(env, share=share)
    assert status == 400
    assert body == {"success": False, "message": "Invalid share ID", "status": "error", "progress": 0}


def test_share_id_cannot_reach_parent_directory(env):
    # "<session>/../status.json" would otherwise read the file directly under the uploads area
    (env.uploads / SESSION).mkdir()
    write_status(env, {"status": "complete", "secret": True}, env.uploads / "status.json")
    status, body = post(env, share="%2E%2E")
    assert status == 400
    assert "secret" not in body


@pytest.mark.parametrize("share", ["a1b2c3d4", "123e4567-e89b-12d3-a456-426614174000", "pytest_1"])
def test_generated_share_ids_are_accepted(env, share):
    write_status(env, {"status": "complete"}, env.uploads / SESSION / share / "status.json")
    status, _ = post(env, share=share)
    assert status == 200


def test_missing_status_file(env):
    status, body = post(env)
    assert status == 404
    assert body == {"success": False, "message": "Job status file not found", "status": "error", "progress": 0}


def test_complete_sets_progress_to_100(env):
    write_status(env, {"status": "complete", "message": "Done", "progress": 40, "dataset_id": "DS1"})
    status, body = post(env)
    assert status == 200
    assert body == {"status": "complete", "message": "Done", "progress": 100, "dataset_id": "DS1"}
    assert env.ps_calls == []


@pytest.mark.parametrize("state", ["uploaded", "extracting", "error", ""])
def test_other_statuses_are_returned_unchanged(env, state):
    write_status(env, {"status": state, "progress": 7, "process_id": 99})
    status, body = post(env)
    assert status == 200
    assert body == {"status": state, "progress": 7, "process_id": 99}
    assert env.ps_calls == []


def test_processing_with_job_id(env):
    write_status(env, {"status": "processing", "progress": 30, "job_id": "J1"})
    status, body = post(env)
    assert status == 200
    assert body == {"status": "processing", "progress": 30, "job_id": "J1"}
    assert env.ps_calls == []


def test_processing_falls_back_to_running_process(env):
    write_status(env, {"status": "processing", "progress": 30, "process_id": 123})
    status, body = post(env)
    assert status == 200
    assert body["status"] == "processing"
    assert body["progress"] == 30
    assert env.ps_calls == [["ps", "-p", "123"]]


def test_processing_with_dead_process(env):
    env.ps.returncode = 1
    write_status(env, {"status": "processing", "progress": 30, "process_id": 123})
    status, body = post(env)
    assert status == 400
    assert body["status"] == "error"
    assert body["progress"] == 0
    assert body["message"] == "The processing step failed. Please contact the gEAR team."
    assert env.ps_calls == [["ps", "-p", "123"]]


def test_processing_without_job_or_process_id(env):
    write_status(env, {"status": "processing", "progress": 30})
    status, body = post(env)
    assert status == 400
    assert body["status"] == "error"
    assert body["progress"] == 0
    assert env.ps_calls == []


def test_gosling_processing_skips_job_checks(env):
    write_status(env, {"status": "processing", "progress": 50, "completed_tracks": 1, "total_tracks": 2})
    status, body = post(env, dataset_format="gosling")
    assert status == 200
    assert body["status"] == "processing"
    assert body["completed_tracks"] == 1
    assert env.ps_calls == []
