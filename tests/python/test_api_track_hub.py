"""POST /import/trackhub/<share_uid>/copy: input validation, auth, and safe handling of uploaded file names."""

import importlib
import importlib.util
import io
import json
import sys
import types
from pathlib import Path

import pytest
from flask import Flask
from flask_restful import Api

FAKES = Path(__file__).resolve().parent / "fakes"
SESSION = "SESS1"
SHARE = "SHARE1"
# www/js/classes/trackhub.js sends each file as "tracks[<track id>][file]"; the resource takes the
#  track id from between the first pair of brackets.
FILE_FIELD = "tracks[t1][file]"


def _load_stub(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def track_hub_module():
    """
    Import resources.track_hub. If gear.trackhub can't be imported here (it needs pyBigWig and
    Biopython), a stub stands in for it while this module's tests run.
    """
    with pytest.MonkeyPatch.context() as mp:
        try:
            importlib.import_module("gear.trackhub")
        except ImportError:
            mp.setitem(sys.modules, "gear.trackhub", _load_stub("gear.trackhub", FAKES / "api_trackhub.py"))
        mp.delitem(sys.modules, "resources.track_hub", raising=False)
        yield importlib.import_module("resources.track_hub")


@pytest.fixture
def env(track_hub_module, monkeypatch, tmp_path):
    """Patch the module for an isolated upload area, one valid session and a disabled queue."""
    mod = track_hub_module
    processed = []

    def fake_process(**kwargs):
        processed.append(kwargs)
        return {"success": True, "message": "processed"}

    monkeypatch.setattr(mod, "user_upload_file_base", tmp_path)
    monkeypatch.setattr(mod.geardb, "get_user_from_session_id",
                        lambda session_id=None: types.SimpleNamespace(id=1) if session_id == SESSION else None)
    monkeypatch.setattr(mod.geardb, "_read_domain_url", lambda: "https://gear.example.org")
    monkeypatch.setattr(mod, "process_trackhub_synchronously", fake_process)
    monkeypatch.delenv("ENVIRONMENT", raising=False)

    # Fresh config with the queue disabled and no higlass section (restored afterwards)
    monkeypatch.setattr(mod, "_config", type(mod._config)())
    mod._config.read_dict({"dataset_uploader": {"queue_enabled": "false"}})

    app = Flask(__name__)
    Api(app).add_resource(mod.TrackHubCopy, "/import/trackhub/<share_uid>/copy")
    client = app.test_client()
    client.set_cookie("gear_session_id", SESSION)
    return types.SimpleNamespace(client=client, tmp=tmp_path, processed=processed, module=mod,
                                 staging=tmp_path / SESSION / SHARE)


def post(env, share=SHARE, hub_json=None, tracks=None, files=None, **extra):
    data = {
        "hub_json": json.dumps({"hub": "myhub", "shortLabel": "My hub"}) if hub_json is None else hub_json,
        "assembly": "mm10",
        "tracks": json.dumps([{"id": "t1", "type": "bigWig", "shortLabel": "Track 1"}]) if tracks is None else tracks,
        "dry_run": "false",
        **extra,
    }
    for field, (filename, content) in (files or {}).items():
        data[field] = (io.BytesIO(content), filename)
    resp = env.client.post(f"/import/trackhub/{share}/copy", data=data, content_type="multipart/form-data")
    return resp.status_code, resp.get_json()


def write_metadata(env, metadata):
    env.staging.mkdir(parents=True, exist_ok=True)
    (env.staging / "metadata.json").write_text(json.dumps(metadata))


def test_invalid_hub_json(env):
    status, body = post(env, hub_json="{not json")
    assert status == 400
    assert "Invalid JSON" in body["message"]


def test_invalid_tracks_json(env):
    status, body = post(env, tracks="[oops")
    assert status == 400
    assert "Invalid JSON" in body["message"]


@pytest.mark.parametrize("share", ["bad share", ".hidden", "bad;share"])
def test_share_id_that_is_not_a_plain_name(env, share):
    status, body = post(env, share=share)
    assert status == 400
    assert body["message"] == "Invalid share ID"
    assert not any(env.tmp.iterdir())


def test_unknown_session(env):
    env.client.set_cookie("gear_session_id", "NOT-A-SESSION")
    status, body = post(env)
    assert status == 401
    assert body["success"] is False
    assert not any(env.tmp.iterdir())


def test_path_traversal_in_file_name_stays_in_staging_area(env):
    write_metadata(env, {"dataset_uid": "D1"})
    post(env, files={FILE_FIELD: ("../../evil.bw", b"EVIL")})
    assert (env.staging / "evil.bw").read_bytes() == b"EVIL"
    found = [p for p in env.tmp.rglob("evil.bw")]
    assert found == [env.staging / "evil.bw"]
    # Nothing escaped tmp_path either (the traversal would have landed next to SESS1/)
    assert not (env.tmp.parent / "evil.bw").exists()
    assert not (env.tmp / "evil.bw").exists()


def test_missing_metadata_returns_json_error_with_job_id(env):
    status, body = post(env, files={FILE_FIELD: ("t1.bw", b"data")})
    assert status == 400
    assert body["success"] is False
    assert "Metadata file not found" in body["message"]
    assert body["job_id"]
    # Nothing was saved, and no staging directory was created for the unknown upload
    assert not env.staging.exists()
    assert not list(env.tmp.rglob("t1.bw"))
    assert env.processed == []


def saved_files(env):
    return sorted(p.name for p in env.staging.iterdir())


@pytest.mark.parametrize("metadata_text, message", [
    (json.dumps({"title": "no uid"}), "Dataset ID not found"),
    ("{not json", "Metadata file could not be read"),
    (json.dumps(["not", "an", "object"]), "Dataset ID not found"),
])
def test_bad_metadata_saves_no_track_files(env, metadata_text, message):
    env.staging.mkdir(parents=True)
    (env.staging / "metadata.json").write_text(metadata_text)
    status, body = post(env, files={FILE_FIELD: ("t1.bw", b"data")})
    assert status == 400
    assert message in body["message"]
    assert body["job_id"]
    assert saved_files(env) == ["metadata.json", "status.json"]
    status_file = json.loads((env.staging / "status.json").read_text())
    assert status_file["status"] == "error" and status_file["job_id"] == body["job_id"]
    assert (env.staging / "metadata.json").read_text() == metadata_text
    assert env.processed == []


def test_missing_domain_url_saves_no_track_files(env, monkeypatch):
    write_metadata(env, {"dataset_uid": "D1"})
    monkeypatch.setattr(env.module.geardb, "_read_domain_url", lambda: "")
    status, body = post(env, files={FILE_FIELD: ("t1.bw", b"data")})
    assert status == 500
    assert "Domain URL not configured" in body["message"]
    assert saved_files(env) == ["metadata.json", "status.json"]
    assert "dataset_format" not in json.loads((env.staging / "metadata.json").read_text())


def test_invalid_file_name_saves_no_track_files(env):
    write_metadata(env, {"dataset_uid": "D1"})
    status, body = post(env, files={"tracks[t0][file]": ("good.bw", b"data"), FILE_FIELD: ("..", b"data")})
    assert status == 400
    assert "Invalid file name" in body["message"]
    assert saved_files(env) == ["metadata.json", "status.json"]


def test_successful_synchronous_import(env):
    write_metadata(env, {"dataset_uid": "D1"})
    status, body = post(env, files={FILE_FIELD: ("my track (1).bw", b"BIGWIG")})
    assert status == 200, body
    assert body["success"] is True
    assert body["message"] == "processed"
    assert body["job_id"]

    assert (env.staging / "my_track_1.bw").read_bytes() == b"BIGWIG"
    assert json.loads((env.staging / "metadata.json").read_text())["dataset_format"] == "gosling"

    (call,) = env.processed
    assert call["job_id"] == body["job_id"]
    assert call["share_uid"] == SHARE
    assert call["staging_area"] == env.staging
    assert call["hub_url"] == "https://gear.example.org/tracks/D1"
    assert call["assembly"] == "mm10"
    assert call["dry_run"] is False
    assert call["higlass_config"] is None
    (stanza,) = call["track_stanzas"]
    assert stanza["id"] == "t1"
    assert stanza["uploadedFileName"] == "my_track_1.bw"
