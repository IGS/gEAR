"""Dataset CGIs: login checks, unknown datasets and single, well-formed JSON responses."""

from pathlib import Path

import pytest

from helpers.cgi_harness import run_cgi

FAKES_DIR = Path(__file__).resolve().parent / "fakes"
ANALYSIS_FAKE = {"gear.analysis": FAKES_DIR / "datasets_analysis.py"}

# Session "owner" is user 1, who owns dataset DS1.
DB = {
    "sessions": {"owner": 1, "other": 2},
    "datasets": [{"id": "DS1", "share_id": "SH1", "owner_id": 1, "title": "A dataset"}],
}


def test_get_dataset_info_unknown_dataset():
    result = run_cgi("get_dataset_info.cgi", db=DB, query={"dataset_id": "NOPE"})
    assert result.json() == {"success": 0, "error": "Dataset not found"}


@pytest.mark.parametrize("script, query, error", [
    ("save_datasetinfo_changes.cgi", {"dataset_id": "DS1", "title": "Hijacked", "visibility": "1"},
     "User must be logged in"),
    ("mark_dataset_for_removal.cgi", {"dataset_id": "DS1"},
     "Not able to remove dataset. User must be logged in."),
    ("save_default_display.cgi", {"dataset_id": "DS1", "display_id": "7", "is_multigene": "0"},
     "User must be logged in"),
])
def test_dataset_changes_require_login(script, query, error):
    result = run_cgi(script, db=DB, query=query)
    body = result.json()
    assert not body["success"] and body["error"] == error
    assert result.writes() == []
    assert not result.logged("dataset.save_change")


class TestValidateShareId:
    # The share ID exists: "SELECT share_id FROM dataset WHERE share_id = %s" echoes it back
    VALID_SHARE = {"match": "FROM dataset WHERE share_id", "rows": "echo"}

    def run(self, sql, **query):
        return run_cgi("validate_share_id.cgi", db={**DB, "sql": sql}, query={"scope": "dataset", **query})

    def test_anonymous_is_asked_to_log_in(self):
        result = self.run([self.VALID_SHARE], share_id="SH1")
        assert result.json() == {"error": "Please log in to add a shared dataset.", "success": 0}
        assert not result.logged("sql", "dataset_shares")

    def test_logged_in_user_reaches_already_shared_check(self):
        result = self.run([self.VALID_SHARE], share_id="SH1", session_id="other")
        assert result.json() == {"success": 1}
        (check,) = result.logged("sql", "dataset_shares")
        assert check["params"] == ["SH1", 2]

    def test_already_shared(self):
        sql = [self.VALID_SHARE, {"match": "JOIN dataset_shares", "rows": [["SH1"]]}]
        body = self.run(sql, share_id="SH1", session_id="other").json()
        assert body["success"] == 0 and "already have this dataset" in body["error"]

    def test_unknown_share_id(self):
        body = self.run([], share_id="NOPE", session_id="other").json()
        assert body["success"] == 0 and "no longer available" in body["error"]


def test_get_dataset_displays_unknown_dataset():
    result = run_cgi("get_dataset_displays.cgi", db=DB, query={"dataset_id": "NOPE", "session_id": "owner"})
    body = result.json()
    assert body["success"] == 0 and body["error"] == "Dataset not found"
    assert body["owner"] == []


def test_copy_dataset_analysis_rejects_crafted_session():
    result = run_cgi("copy_dataset_analysis.cgi", db=DB, fake_modules=ANALYSIS_FAKE, query={
        "session_id": "../..", "dataset_id": "DS1",
        "source_analysis_id": "A1", "source_analysis_type": "public",
        "dest_analysis_id": "A2", "dest_analysis_type": "user_unsaved",
    })
    assert result.json() == {"success": 0, "error": "You must be logged in to copy an analysis."}
    assert "FAKE get_analysis" not in result.stderr and "FAKE Analysis" not in result.stderr


class TestCreateAccount:
    FORM = {"first-last": "Some One", "email": "someone@example.org", "password": "pw",
            "verification_code_long": "anything-long", "verification_code_short": "SHORT"}

    def test_existing_email_gets_one_error(self):
        db = {"sql": [{"match": "FROM guser", "rows": [[5]]}]}
        result = run_cgi("create_account.cgi", db=db, query=self.FORM)
        body = result.json()    # fails if a second, "success" object was printed
        assert body["success"] == 0 and body["error"] == "User already exists"
        assert body["session_id"] == -1
        assert result.writes() == []

    def test_verification_mismatch(self):
        result = run_cgi("create_account.cgi", query={**self.FORM, "verification_code_short": "WRONG"})
        body = result.json()
        assert body["success"] == 0 and "Verification code mismatch" in body["error"]
        assert result.writes() == []

    def test_new_user_is_created(self):
        result = run_cgi("create_account.cgi", query=self.FORM)
        body = result.json()
        assert body["success"] == 1 and body["session_id"] not in (0, -1)
        assert [w["query"].split()[2] for w in result.writes()] == ["guser", "user_session"]
        assert result.logged("commit")


class TestSaveDatasetDisplay:
    """Saving a display names its preview image by the new display's ID (issue #490)."""

    FAKE_REQUESTS = {"requests": FAKES_DIR / "display_requests.py"}
    CONFIG = '{"gene_symbol": "Sox2"}'

    def run(self, db=None, **query):
        query = {"session_id": "owner", "dataset_id": "DS1", "plot_type": "bar", "plotly_config": self.CONFIG, **query}
        result = run_cgi("save_dataset_display.cgi", query=query, db={**DB, **(db or {})},
                         fake_modules=self.FAKE_REQUESTS)
        assert result.returncode == 0, result.stderr[-2000:]
        return result

    @pytest.mark.parametrize("label", ["", "My display"])
    def test_new_display_uses_inserted_id(self, label):
        result = self.run(db={"lastrowid": 42}, label=label)
        assert result.json() == {"display_id": 42, "success": True}
        (insert,) = result.writes()
        assert insert["query"].startswith("INSERT INTO dataset_display")
        # The preview is requested for this dataset and named after display 42, never "None"
        assert [e["url"] for e in result.logged("requests.post")][0] == "https://localhost/api/plot/DS1"
        assert "display id 42" in result.stderr
        assert "None" not in result.stderr

    def test_update_by_owner_regenerates_preview(self):
        result = self.run(db={"sql": [{"match": "SELECT user_id FROM dataset_display", "rows": [[1]]}]},
                          id="7", label="Renamed")
        assert result.json() == {"display_id": "7", "success": True}
        assert [w["query"].split()[0] for w in result.writes()] == ["UPDATE"]
        assert result.logged("requests.post")

    def test_update_by_other_user_changes_nothing(self):
        result = self.run(db={"sql": [{"match": "SELECT user_id FROM dataset_display", "rows": [[2]]}]},
                          id="7", label="Hijacked")
        assert result.json() == {"display_id": "7", "success": False}
        assert result.writes() == []
        # The owner's preview image is not overwritten with this user's config
        assert not result.logged("requests.post")

    def test_unknown_display_id(self):
        result = self.run(id="999", label="x")
        assert result.json() == {"display_id": "999", "success": False}
        assert not result.logged("requests.post")


class TestSaveDefaultDisplay:
    def run(self, **query):
        query = {"session_id": "owner", "dataset_id": "DS1", "is_multigene": "0", **query}
        return run_cgi("save_default_display.cgi", query=query, db=DB)

    @pytest.mark.parametrize("display_id", [None, "", "None", "undefined", "7; DROP"])
    def test_invalid_display_id(self, display_id):
        result = self.run(**({} if display_id is None else {"display_id": display_id}))
        assert result.json() == {"success": False, "error": "Invalid display ID"}
        assert result.writes() == []

    def test_valid_display_id_saves_preference(self):
        result = self.run(display_id="7")
        assert result.json() == {"success": True}
        (insert,) = result.writes()
        assert "INSERT INTO dataset_preference" in insert["query"]
        assert insert["params"] == [1, "DS1", "7", 0]

    def test_failed_insert_skips_symlink(self):
        result = run_cgi("save_default_display.cgi", db={**DB, "fail_sql": ["INSERT INTO dataset_preference"]},
                         query={"session_id": "owner", "dataset_id": "DS1", "display_id": "7", "is_multigene": "0"})
        assert result.returncode == 0, result.stderr[-2000:]
        assert result.json()["success"] is False
        assert not result.logged("sql", "dataset_preference dp")
