"""Account-related CGIs (preferences, email check, user history, notes): login checks and one JSON response."""

from pathlib import Path

import pytest

from helpers.cgi_harness import run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"
USERHISTORY = {"gear.userhistory": FAKES / "accounts_userhistory.py"}

# Session "owner" is user 1; "other" is user 2.
SESSIONS = {"sessions": {"owner": 1, "other": 2}}


class TestSaveUserDefaultOrganism:
    def test_anonymous_is_refused(self):
        result = run_cgi("save_user_default_organism.cgi", db=SESSIONS, query={"default_org_id": "1"})
        assert result.json() == {"success": False, "error": "User not logged in"}
        assert result.writes() == []

    def test_logged_in_user_saves(self):
        result = run_cgi("save_user_default_organism.cgi", db=SESSIONS,
                         query={"session_id": "owner", "default_org_id": "3"})
        assert result.json() == {"success": True}
        (update,) = result.writes()
        assert update["query"].startswith("UPDATE guser SET default_org_id")
        assert update["params"] == ["3", 1]


class TestCheckExistingEmail:
    def test_no_email(self):
        body = run_cgi("check_existing_email.cgi").json()
        assert body["email_exists"] == 0
        assert body["error"]

    def test_existing_email(self):
        result = run_cgi("check_existing_email.cgi", query={"email": " someone@example.org "},
                         db={"sql": [{"match": "FROM guser", "rows": [[1]]}]})
        assert result.json() == {"email_exists": 1}
        (query,) = result.logged("sql", "FROM guser")
        assert query["params"] == ["someone@example.org"]     # whitespace removed

    def test_unknown_email(self):
        result = run_cgi("check_existing_email.cgi", query={"email": "nobody@example.org"})
        assert result.json() == {"email_exists": 0}


class TestAddToUserHistory:
    def test_anonymous_is_refused(self):
        result = run_cgi("add_to_user_history.cgi", db=SESSIONS, fake_modules=USERHISTORY,
                         query={"entry_category": "dataset_search", "label": "x"})
        assert result.json() == {"success": 0, "error": "User must be logged in"}
        assert not result.logged("userhistory.add_record")

    def test_logged_in_user_adds_record(self):
        result = run_cgi("add_to_user_history.cgi", db=SESSIONS, fake_modules=USERHISTORY,
                         query={"session_id": "owner", "entry_category": "dataset_search", "label": "Sox2"})
        assert result.json() == {"success": 1}
        (record,) = result.logged("userhistory.add_record")
        assert record["user_id"] == 1
        assert record["entry_category"] == "dataset_search"
        assert record["label"] == "Sox2"


class TestApplyNoteChanges:
    # Note 5 belongs to user 1
    DB = {**SESSIONS, "sql": [{"match": "FROM note", "rows": [[1]]}]}
    EDIT = {"scope": "edit", "note_id": "5", "title": "T", "ldesc": "D", "access_level": "1"}

    def test_owner_can_edit(self):
        result = run_cgi("apply_note_changes.cgi", db=self.DB, query={"session_id": "owner", **self.EDIT})
        assert result.json()["success"] == 1
        (update,) = result.writes()
        assert update["query"].startswith("UPDATE note")
        assert update["params"] == ["T", "D", "1", "5"]

    def test_other_user_cannot_edit(self):
        result = run_cgi("apply_note_changes.cgi", db=self.DB, query={"session_id": "other", **self.EDIT})
        body = result.json()
        assert body["success"] == 0
        assert "Not note owner." in body["error"]
        assert not [w for w in result.writes() if w["query"].startswith("UPDATE")]

    @pytest.mark.parametrize("scope", ["edit", "remove", "new"])
    def test_anonymous_is_refused(self, scope):
        result = run_cgi("apply_note_changes.cgi", db=self.DB, query={**self.EDIT, "scope": scope})
        body = result.json()
        assert body["success"] == 0 and "log in" in body["error"]
        assert result.writes() == []

    def test_edit_without_note_id(self):
        query = {k: v for k, v in self.EDIT.items() if k != "note_id"}
        result = run_cgi("apply_note_changes.cgi", db=self.DB, query={"session_id": "owner", **query})
        assert "Invalid note ID" in result.json()["error"]
        assert result.writes() == []

    def test_other_user_cannot_remove(self):
        result = run_cgi("apply_note_changes.cgi", db=self.DB,
                         query={"session_id": "other", "scope": "remove", "note_id": "5"})
        assert "Not note owner." in result.json()["error"]
        assert result.writes() == []
