"""Password hashing in gear.passwords (issue #1350)."""

import base64
import hashlib
import re
from pathlib import Path

import pytest
from gear import passwords

REPO_ROOT = Path(__file__).resolve().parents[2]

PASSWORD = "correct horse"
SCRYPT = passwords.hash_password(PASSWORD)
SHA3 = hashlib.sha3_256(PASSWORD.encode()).hexdigest()
MD5 = hashlib.md5(PASSWORD.encode()).hexdigest()


def make_hash(n=passwords.SCRYPT_N, r=passwords.SCRYPT_R, p=passwords.SCRYPT_P, salt=b"s" * 16, key=b"k" * 32):
    return "$".join(["scrypt", str(n), str(r), str(p), base64.b64encode(salt).decode(), base64.b64encode(key).decode()])


class TestHashPassword:
    def test_format_records_the_cost_parameters(self):
        algorithm, n, r, p, salt, key = SCRYPT.split("$")
        assert (algorithm, int(n), int(r), int(p)) == ("scrypt", passwords.SCRYPT_N, passwords.SCRYPT_R, passwords.SCRYPT_P)
        assert len(base64.b64decode(salt)) == passwords.SALT_BYTES
        assert len(base64.b64decode(key)) == passwords.KEY_BYTES

    def test_same_password_gets_a_different_hash_each_time(self):
        """Salted, so two users with the same password (or one user setting it twice) don't share a hash."""
        assert passwords.hash_password(PASSWORD) != SCRYPT

    def test_matches_the_scrypt_reference_implementation(self):
        _, n, r, p, salt, key = SCRYPT.split("$")
        expected = hashlib.scrypt(PASSWORD.encode(), salt=base64.b64decode(salt), n=int(n), r=int(r), p=int(p),
                                  dklen=32, maxmem=2**28)
        assert base64.b64decode(key) == expected

    @pytest.mark.parametrize("password", ["", "pässwörd ☃", "x" * 5000, "with$dollar$signs"])
    def test_round_trip(self, password):
        encoded = passwords.hash_password(password)
        assert passwords.verify_password(password, encoded)
        assert not passwords.verify_password(password + "x", encoded)

    def test_fits_in_the_database_column(self):
        schema = (REPO_ROOT / "create_schema.sql").read_text()
        column_length = int(re.search(r"^\s*pass\s+VARCHAR\((\d+)\)", schema, re.MULTILINE).group(1))
        assert len(SCRYPT) <= passwords.MAX_HASH_LENGTH <= column_length


class TestVerifyPassword:
    @pytest.mark.parametrize("stored", [SCRYPT, SHA3, MD5, SHA3.upper()], ids=["scrypt", "sha3", "md5", "sha3-upper"])
    def test_correct_password(self, stored):
        assert passwords.verify_password(PASSWORD, stored)

    @pytest.mark.parametrize("stored", [SCRYPT, SHA3, MD5], ids=["scrypt", "sha3", "md5"])
    @pytest.mark.parametrize("password", ["wrong", "", None, PASSWORD + " ", PASSWORD.upper(), SCRYPT, SHA3, MD5])
    def test_wrong_password(self, password, stored):
        assert not passwords.verify_password(password, stored)

    @pytest.mark.parametrize("stored", [
        None, "", "plaintext", PASSWORD, SHA3[:-1], MD5 + "0", "z" * 64,
        SCRYPT[:64],                                  # truncated by a VARCHAR(64) column
        SCRYPT.replace("scrypt$", "bcrypt$"),
        SCRYPT + "$extra",
        make_hash(n=0), make_hash(n=1), make_hash(n=1000), make_hash(n=-(2**14)), make_hash(n="x"),
        make_hash(n=passwords.MAX_SCRYPT_N * 2),      # would need over 2 GB
        make_hash(r=0), make_hash(r=passwords.MAX_SCRYPT_R + 1),
        make_hash(p=0), make_hash(p=passwords.MAX_SCRYPT_P + 1),
        make_hash(salt=b""), make_hash(key=b""),
        "scrypt$16384$8$1$not base64!$also not",
    ])
    def test_unusable_stored_value_never_matches(self, stored):
        for password in (PASSWORD, "", stored):
            assert not passwords.verify_password(password, stored)

    def test_hash_made_with_other_settings_still_verifies(self, monkeypatch):
        """Raising the cost later must not invalidate the hashes already stored."""
        monkeypatch.setattr(passwords, "SCRYPT_N", 2**14)
        weaker = passwords.hash_password(PASSWORD)
        monkeypatch.undo()
        assert weaker.startswith("scrypt$16384$")
        assert passwords.verify_password(PASSWORD, weaker)
        assert not passwords.verify_password("wrong", weaker)


class TestNeedsRehash:
    def test_current_hash_is_kept(self):
        assert not passwords.needs_rehash(SCRYPT)

    @pytest.mark.parametrize("stored", [
        SHA3, MD5, None, "",
        make_hash(n=2**14), make_hash(r=4), make_hash(p=2), make_hash(salt=b"s" * 8), make_hash(key=b"k" * 16),
    ])
    def test_older_or_weaker_hash_is_replaced(self, stored):
        assert passwords.needs_rehash(stored)


def test_no_script_hashes_passwords_directly():
    """Everything that stores or checks a password goes through gear.passwords."""
    offenders = [path.name for path in sorted((REPO_ROOT / "www" / "cgi").glob("*.cgi"))
                 if re.search(r"hashlib\.(md5|sha\w+)\(\s*\w*pass", path.read_text(errors="replace"))]
    assert offenders == []
