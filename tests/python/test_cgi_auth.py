"""Login and password reset CGIs, plus the e-mail verification code helpers in geardb."""

import hashlib
import uuid

import pytest
from gear import passwords
from helpers.cgi_harness import run_cgi

PASSWORD = "correct horse"
SHA3 = hashlib.sha3_256(PASSWORD.encode()).hexdigest()
MD5 = hashlib.md5(PASSWORD.encode()).hexdigest()
SCRYPT = passwords.hash_password(PASSWORD)

# What the database reports once a new hash has been stored in full (see upgrade_stored_password in login.v2.cgi)
STORED_IN_FULL = {"match": "SELECT CHAR_LENGTH(pass) FROM guser", "rows": [[len(SCRYPT)]]}


class TestUpdatePassword:
    def run(self, count=1, **form):
        db = {"sql": [{"match": "SELECT COUNT(*) FROM guser", "rows": [[count]]}]}
        return run_cgi("update_password.cgi", form=form, db=db)

    def test_valid_help_id_updates_password_and_replaces_help_id(self):
        """The help_id is in the reset e-mail, so it must change once used or the link works forever."""
        result = self.run(password="n3w-pass", help_id="abc123")
        assert result.json() == {"error": "", "success": 1}
        (update,) = result.writes()
        assert update["query"] == "UPDATE guser SET pass = %s, help_id = %s WHERE help_id = %s"
        encoded_pass, new_help_id, old_help_id = update["params"]
        assert encoded_pass.startswith("scrypt$") and passwords.verify_password("n3w-pass", encoded_pass)
        assert not passwords.needs_rehash(encoded_pass)
        assert old_help_id == "abc123"
        assert new_help_id != "abc123" and str(uuid.UUID(new_help_id)) == new_help_id
        assert result.logged("commit")

    def test_each_reset_gets_a_different_help_id(self):
        first = self.run(password="n3w-pass", help_id="abc123").writes()[0]["params"][1]
        second = self.run(password="n3w-pass", help_id="abc123").writes()[0]["params"][1]
        assert first != second

    @pytest.mark.parametrize("form", [{}, {"help_id": ""}])
    def test_missing_help_id(self, form):
        result = self.run(password="n3w-pass", **form)
        assert result.json() == {"error": "No help ID provided.", "success": 0}
        assert result.writes() == []

    @pytest.mark.parametrize("password", [None, ""])
    def test_missing_or_empty_password(self, password):
        form = {"help_id": "abc123"}
        if password is not None:
            form["password"] = password
        result = self.run(**form)
        assert result.json() == {"error": "No password provided.", "success": 0}
        assert result.writes() == []

    def test_unknown_help_id(self):
        result = self.run(count=0, password="n3w-pass", help_id="nobody")
        assert result.json() == {"error": "Help ID not found.", "success": 0}
        assert result.writes() == []

    def test_ambiguous_help_id(self):
        result = self.run(count=2, password="n3w-pass", help_id="abc")
        assert result.json() == {"error": "Duplicate help ID prefixes", "success": 0}
        assert result.writes() == []


class TestLogin:
    def run(self, stored_hash=SCRYPT, password=PASSWORD, extra_sql=(), user_found=True):
        rows = [[7, "Alice", stored_hash, 0]] if user_found else []
        db = {"sql": [{"match": "FROM guser WHERE email", "rows": rows}, *extra_sql]}
        form = {"user-email": "alice@example.org"}
        if password is not None:
            form["user-password"] = password
        return run_cgi("login.v2.cgi", form=form, db=db)

    def assert_logged_in(self, result):
        body = result.json()
        (insert,) = [w for w in result.writes() if "guser" not in w["query"]]
        assert insert["query"] == "INSERT INTO user_session (user_id, session_id) VALUES (%s, %s)"
        assert insert["params"][0] == 7
        assert body == {"session_id": insert["params"][1], "is_admin": 0, "user_name": "Alice",
                        "gear_default_domain": "Hearing (site default)", "default_profile_share_id": 0}

    def test_scrypt_password_logs_in_and_is_left_alone(self):
        result = self.run()
        self.assert_logged_in(result)
        assert not [w for w in result.writes() if "guser" in w["query"]]

    @pytest.mark.parametrize("stored_hash", [SHA3, MD5], ids=["sha3", "md5"])
    def test_legacy_hash_logs_in_and_is_upgraded(self, stored_hash):
        """Issue #1350: unsalted hashes are replaced at login, the only time the plaintext is available."""
        result = self.run(stored_hash=stored_hash, extra_sql=[STORED_IN_FULL])
        self.assert_logged_in(result)
        (update,) = [w for w in result.writes() if "guser" in w["query"]]
        assert update["query"] == "UPDATE guser SET pass = %s WHERE id = %s AND pass = %s"
        new_hash, user_id, old_hash = update["params"]
        assert (user_id, old_hash) == (7, stored_hash)
        assert passwords.verify_password(PASSWORD, new_hash) and not passwords.needs_rehash(new_hash)
        assert not result.logged("rollback")

    def test_weaker_scrypt_settings_are_upgraded(self, monkeypatch):
        monkeypatch.setattr(passwords, "SCRYPT_N", 2**14)
        weak = passwords.hash_password(PASSWORD)
        monkeypatch.undo()
        result = self.run(stored_hash=weak, extra_sql=[STORED_IN_FULL])
        self.assert_logged_in(result)
        (update,) = [w for w in result.writes() if "guser" in w["query"]]
        assert not passwords.needs_rehash(update["params"][0])

    @pytest.mark.parametrize("stored_length", [64, None], ids=["truncated", "row-missing"])
    def test_upgrade_is_rolled_back_if_the_new_hash_was_not_stored_in_full(self, stored_length):
        """
        If guser.pass is still VARCHAR(64) and MySQL isn't strict, the new hash is cut short. Keeping
        it would lock the user out, so it is rolled back and the login still succeeds.
        """
        rows = [[stored_length]] if stored_length else []
        result = self.run(stored_hash=SHA3, extra_sql=[{"match": "SELECT CHAR_LENGTH(pass) FROM guser", "rows": rows}])
        self.assert_logged_in(result)
        assert len(result.logged("rollback")) == 1
        assert "could not upgrade the password hash of user 7" in result.stderr
        # The rollback comes before the session is added, so the session is still committed
        events = [e["event"] if e["event"] != "sql" else e["query"].split()[0] for e in result.db_log]
        assert events.index("rollback") < events.index("INSERT") < len(events) - 1 - events[::-1].index("commit")

    def test_upgrade_failure_does_not_block_login(self):
        """MySQL in strict mode raises instead of truncating."""
        strict = run_cgi("login.v2.cgi", form={"user-email": "alice@example.org", "user-password": PASSWORD},
                         db={"sql": [{"match": "FROM guser WHERE email", "rows": [[7, "Alice", MD5, 0]]}],
                             "fail_sql": ["UPDATE guser SET pass"]})
        self.assert_logged_in(strict)
        assert len(strict.logged("rollback")) == 1

    def test_existing_session_and_profile_are_reused(self):
        result = self.run(extra_sql=[{"match": "FROM user_session", "rows": [["existing-session"]]},
                                     {"match": "FROM layout", "rows": [["My profile", "LAYOUT1"]]}])
        body = result.json()
        assert body["session_id"] == "existing-session"
        assert body["gear_default_domain"] == "My profile" and body["default_profile_share_id"] == "LAYOUT1"
        assert result.writes() == []

    @pytest.mark.parametrize("stored_hash", [SCRYPT, SHA3, MD5], ids=["scrypt", "sha3", "md5"])
    @pytest.mark.parametrize("password", ["wrong", "", None, SHA3, MD5, SCRYPT])
    def test_wrong_password(self, password, stored_hash):
        """Including the stored hash itself, which must not work as a password."""
        result = self.run(password=password, stored_hash=stored_hash)
        assert result.json() == {"session_id": -1, "is_admin": 0}
        assert result.writes() == []

    def test_unknown_user(self):
        result = self.run(user_found=False)
        assert result.json() == {"session_id": -1, "is_admin": 0}
        assert result.writes() == []
        (query,) = result.logged("sql")
        assert query["params"] == ["alice@example.org"]

    def test_unknown_user_and_wrong_password_look_the_same(self):
        """The response must not reveal whether an e-mail is registered; the reason only goes to the server log."""
        unknown = self.run(user_found=False)
        wrong = self.run(password="wrong")
        assert unknown.json() == wrong.json()
        assert unknown.writes() == wrong.writes() == []
        assert "e-mail not found" in unknown.stderr
        assert "incorrect password" in wrong.stderr


class TestVerificationCode:
    def test_short_form_is_first_character_of_each_part(self):
        import geardb
        assert geardb.get_verification_code_short_form("a1b2-c3d4-e5f6-0789") == "ace0"
        assert geardb.get_verification_code_short_form("nodash") == "n"

    @pytest.mark.parametrize("short_form, expected", [("ace0", True), ("ACE0", False), ("ace", False), ("", False)])
    def test_check(self, short_form, expected):
        import geardb
        assert geardb.check_verification_code("a1b2-c3d4-e5f6-0789", short_form) is expected
