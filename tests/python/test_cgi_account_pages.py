"""Account page CGIs: session info, profile changes, password-reset links and the user's history."""

import hashlib
import uuid
from pathlib import Path

import pytest
from helpers.cgi_harness import run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"
USERHISTORY = {"gear.userhistory": FAKES / "cgi_acct_userhistory.py"}

# Session "owner" is user 1; "other" is user 2.
SESSIONS = {"owner": 1, "other": 2}

# email, user_name, is_admin, id, institution, colorblind_mode, updates_wanted (TestGetSessionInfoV2 adds
#  default_org_id and the layout share_id)
GUSER_ROW = ["sally@example.org", "Sally", 0, 1, "UMB", 1, 0]


def assert_clean_exit(result):
    assert result.returncode == 0, result.stderr[-2000:]
    assert "Traceback" not in result.stderr, result.stderr[-2000:]


# ---------------------------------------------------------------- get_session_info.v2.cgi

class TestGetSessionInfoV2:
    def run(self, layout_share_id="L1", **db):
        row = [*GUSER_ROW, 3, layout_share_id]
        return run_cgi("get_session_info.v2.cgi", query={"session_id": "owner"},
                       db={"sql": [{"match": "FROM guser u", "rows": [row]}], **db})

    def test_known_session(self):
        result = self.run()
        assert result.json() == {"email": "sally@example.org", "user_name": "Sally", "is_admin": 0,
                                 "institution": "UMB", "colorblind_mode": 1, "updates_wanted": 0,
                                 "default_org_id": 3, "layout_share_id": "L1", "success": 1}
        assert result.writes() == []

    def test_no_saved_layout_uses_site_default(self):
        result = self.run(layout_share_id=None, servercfg={"content": {"default_layout_share_id": "SITE"}})
        assert result.json()["layout_share_id"] == "SITE"

    @pytest.mark.parametrize("session", ["expired", None])
    def test_unknown_session(self, session):
        query = {"session_id": session} if session else {}
        result = run_cgi("get_session_info.v2.cgi", query=query)
        assert result.json() == {"email": None, "user_name": None, "success": 0}


# ---------------------------------------------------------------- save_user_account_changes.cgi

class TestSaveUserAccountChanges:
    DB = {"sessions": SESSIONS}

    def run(self, **form):
        return run_cgi("save_user_account_changes.cgi", db=self.DB, form=form)

    @pytest.mark.parametrize("session, error", [(None, "No session id"), ("expired", "Invalid session id")])
    @pytest.mark.parametrize("scope", ["password", "settings_change"])
    def test_anonymous_is_refused(self, session, error, scope):
        form = {"scope": scope, "new_password": "hunter22", "help_id": "H1", "email": "x@example.org"}
        if session:
            form["session_id"] = session
        result = self.run(**form)
        body = result.json()
        assert body["success"] == 0 and error in body["error"]
        assert result.writes() == []

    def test_password_change_without_password(self):
        result = self.run(session_id="owner", scope="password", help_id="H1")
        assert result.json() == {"error": "No new password provided.", "success": 0}
        assert result.writes() == []

    def test_password_change_rotates_help_id(self):
        """The password-reset help_id is replaced after the password is changed, so an old reset link stops working."""
        result = self.run(session_id="owner", scope="password", help_id="H1", new_password="hunter22")
        assert result.json() == {"success": 1}
        password, help_id = result.writes()
        assert password["query"] == "UPDATE guser SET pass = %s WHERE help_id = %s"
        assert password["params"] == [hashlib.sha3_256(b"hunter22").hexdigest(), "H1"]
        assert help_id["query"] == "UPDATE guser SET help_id = %s WHERE help_id = %s"
        new_help_id, old_help_id = help_id["params"]
        assert old_help_id == "H1"
        assert str(uuid.UUID(new_help_id)) == new_help_id
        assert len(result.logged("commit")) == 2

    def test_settings_change(self):
        result = self.run(session_id="owner", scope="settings_change", email="new@example.org", institution="JHU",
                          colorblind_mode="true", want_updates="false", new_password="hunter22")
        assert result.json() == {"success": 1}
        (update,) = result.writes()
        assert update["query"] == ("UPDATE guser SET email = %s, institution = %s, colorblind_mode = %s, "
                                   "updates_wanted = %s, pass = %s, help_id = %s WHERE id = %s")
        *params, new_help_id, user_id = update["params"]
        assert params == ["new@example.org", "JHU", 1, 0, hashlib.sha3_256(b"hunter22").hexdigest()]
        assert user_id == 1 and str(uuid.UUID(new_help_id)) == new_help_id

    def test_settings_change_keeps_password_when_blank(self):
        result = self.run(session_id="other", scope="settings_change", colorblind_mode="false", want_updates="true",
                          new_password="")
        assert result.json() == {"success": 1}
        (update,) = result.writes()
        assert update["query"] == "UPDATE guser SET colorblind_mode = %s, updates_wanted = %s WHERE id = %s"
        assert update["params"] == [0, 1, 2]

    def test_settings_password_change_rotates_help_id(self):
        """Changing the password from the profile page also replaces help_id, so an earlier reset link stops working."""
        result = self.run(session_id="owner", scope="settings_change", new_password="hunter22")
        assert result.json() == {"success": 1}
        (update,) = result.writes()
        assert "pass = %s" in update["query"] and "help_id = %s" in update["query"]
        new_help_id = update["params"][-2]
        assert str(uuid.UUID(new_help_id)) == new_help_id

    def test_email_change_checks_other_accounts(self):
        """The new email address is looked up case-insensitively among the other accounts, and stored trimmed."""
        result = self.run(session_id="owner", scope="settings_change", email="  New@Example.org ")
        assert result.json() == {"success": 1}
        (lookup,) = result.logged("sql", "FROM guser WHERE")
        assert lookup["query"] == "SELECT id FROM guser WHERE LOWER(email) = LOWER(%s) AND id != %s"
        assert lookup["params"] == ["New@Example.org", 1]
        (update,) = result.writes()
        assert update["params"][0] == "New@Example.org"

    def test_email_used_by_another_account_is_refused(self):
        result = run_cgi("save_user_account_changes.cgi", form={
            "session_id": "owner", "scope": "settings_change", "email": "TAKEN@example.org", "institution": "JHU",
            "new_password": "hunter22"}, db={**self.DB, "sql": [{"match": "FROM guser WHERE LOWER(email)",
                                                                  "rows": [[2]]}]})
        body = result.json()
        assert body["success"] == 0 and "already uses this email" in body["error"]
        assert result.writes() == []
        assert not result.logged("commit")

    @pytest.mark.parametrize("email", ["not-an-email", "sally@", "@example.org", "sally@example", "sal ly@example.org",
                                       "sally@example.org, x@example.org", "sally@exa_mple.org"])
    def test_invalid_email_is_refused(self, email):
        result = self.run(session_id="owner", scope="settings_change", email=email, colorblind_mode="true")
        body = result.json()
        assert body["success"] == 0 and "valid email" in body["error"]
        assert result.writes() == []
        assert not result.logged("sql", "FROM guser WHERE")

    def test_current_email_is_not_rechecked(self):
        """The profile page always sends the email, so an unchanged one (even one that predates the format
        check) doesn't block saving other settings."""
        db = {**self.DB, "users": {"1": {"email": "Sally@Lab"}}}
        result = run_cgi("save_user_account_changes.cgi", db=db,
                         form={"session_id": "owner", "scope": "settings_change", "email": "sally@lab ",
                               "institution": "JHU"})
        assert result.json() == {"success": 1}
        assert not result.logged("sql", "FROM guser WHERE")
        (update,) = result.writes()
        assert update["params"][:2] == ["sally@lab", "JHU"]

    @pytest.mark.parametrize("scope", [None, "email"])
    def test_unknown_scope(self, scope):
        form = {"session_id": "owner"}
        if scope:
            form["scope"] = scope
        result = self.run(**form)
        assert result.json()["success"] == 0
        assert result.writes() == []


# ---------------------------------------------------------------- validate_help_id.cgi

class TestValidateHelpId:
    HELP_ID = "3f2a9c1e-7b4d-4e1a-9c2b-5d6e7f809a1b"

    def test_valid_help_id(self):
        result = run_cgi("validate_help_id.cgi", query={"help_id": self.HELP_ID},
                         db={"sql": [{"match": "FROM guser", "rows": [["Sally"]]}]})
        assert result.json() == {"user_name": "Sally", "success": 1}
        (query,) = result.logged("sql", "FROM guser")
        assert query["params"] == [self.HELP_ID]
        assert result.writes() == []

    def test_unknown_help_id(self):
        result = run_cgi("validate_help_id.cgi", query={"help_id": self.HELP_ID})
        body = result.json()
        assert body["success"] == 0 and "expired" in body["error"]
        assert result.writes() == []

    @pytest.mark.parametrize("known", [True, False])
    def test_help_id_is_not_logged(self, known):
        """The help_id is a password-reset token, so it must not be written to the server's error log."""
        db = {"sql": [{"match": "FROM guser", "rows": [["Sally"]]}]} if known else {}
        result = run_cgi("validate_help_id.cgi", query={"help_id": self.HELP_ID}, db=db)
        assert_clean_exit(result)
        assert self.HELP_ID not in result.stderr


# ---------------------------------------------------------------- get_user_history_entries.cgi

class TestGetUserHistoryEntries:
    ROWS = [["2026-10-01 10:00:00", "gene_search", "Sox2", "/p?g=Sox2"],
            ["2026-09-30 09:00:00", "dataset_search", "utricle", "/d?s=utricle"]]
    DB = {"sessions": SESSIONS, "sql": [{"match": "FROM user_history", "rows": ROWS}]}

    def run(self, **query):
        return run_cgi("get_user_history_entries.cgi", db=self.DB, fake_modules=USERHISTORY, query=query)

    def history_query(self, result):
        (entry,) = result.logged("sql", "FROM user_history")
        return entry

    def test_anonymous_gets_nothing(self):
        result = self.run()
        assert result.json() == []
        assert not result.logged("sql", "FROM user_history")

    def test_logged_in_user_gets_latest_entries(self):
        result = self.run(session_id="owner")
        body = result.json()
        assert [(e["entry_category"], e["label"], e["url"]) for e in body] == [
            ("gene_search", "Sox2", "/p?g=Sox2"), ("dataset_search", "utricle", "/d?s=utricle")]
        entry = self.history_query(result)
        assert entry["query"].endswith("ORDER BY entry_date DESC LIMIT %s")
        assert entry["params"] == [1, 5]

    def test_num_entries(self):
        result = self.run(session_id="owner", num_entries="12")
        assert self.history_query(result)["params"] == [1, 12]

    @pytest.mark.parametrize("num_entries", ["5; DELETE FROM guser", "5 UNION SELECT pass FROM guser", "x", "-1", "0"])
    def test_num_entries_must_be_a_positive_number(self, num_entries):
        """num_entries goes into the LIMIT clause, so anything but a positive number falls back to the default."""
        result = self.run(session_id="owner", num_entries=num_entries)
        assert_clean_exit(result)
        assert isinstance(result.json(), list)
        assert self.history_query(result)["params"] == [1, 5]
