"""
Password hashing for user accounts (the `guser.pass` column).

New passwords are hashed with scrypt, which is salted, deliberately slow and memory-hard, using
only the standard library (`hashlib.scrypt`). The stored value records the algorithm and its cost
parameters, so they can be raised later without invalidating existing passwords:

    scrypt$<n>$<r>$<p>$<salt, base64>$<derived key, base64>

Accounts created before this was introduced store an unsalted SHA3-256 or MD5 hexdigest. Those are
still accepted by `verify_password()`; `needs_rehash()` reports them so that login can replace them
with a scrypt hash while it has the plaintext password (see www/cgi/login.v2.cgi).
"""

import base64
import hashlib
import hmac
import os

ALGORITHM = "scrypt"

# OWASP's recommended minimum for scrypt. Memory use is about 128 * n * r bytes (128 MB here).
SCRYPT_N = 2**17
SCRYPT_R = 8
SCRYPT_P = 1
SALT_BYTES = 16
KEY_BYTES = 32

# Upper bounds for parameters read back from the database, so a corrupted or tampered row can't
#  make a login attempt use an unreasonable amount of memory or CPU
MAX_SCRYPT_N = 2**20
MAX_SCRYPT_R = 16
MAX_SCRYPT_P = 4

# guser.pass must hold at least this many characters (see create_schema.sql)
MAX_HASH_LENGTH = 255

LEGACY_HEX_LENGTHS = {
    64: hashlib.sha3_256,
    32: hashlib.md5,
}


def _scrypt(password: str, salt: bytes, n: int, r: int, p: int, dklen: int) -> bytes:
    # OpenSSL refuses to use more than 32 MB unless maxmem is raised
    maxmem = 128 * n * r * 2
    return hashlib.scrypt(password.encode("utf-8"), salt=salt, n=n, r=r, p=p, dklen=dklen, maxmem=maxmem)


def hash_password(password: str) -> str:
    """Return the value to store in `guser.pass` for a new password."""
    salt = os.urandom(SALT_BYTES)
    key = _scrypt(password, salt, SCRYPT_N, SCRYPT_R, SCRYPT_P, KEY_BYTES)
    return "$".join([
        ALGORITHM, str(SCRYPT_N), str(SCRYPT_R), str(SCRYPT_P),
        base64.b64encode(salt).decode("ascii"), base64.b64encode(key).decode("ascii"),
    ])


def _parse(stored_hash: str):
    """Return (n, r, p, salt, key) for a scrypt hash made by `hash_password()`, or None."""
    parts = stored_hash.split("$")
    if len(parts) != 6 or parts[0] != ALGORITHM:
        return None
    try:
        n, r, p = (int(part) for part in parts[1:4])
        salt = base64.b64decode(parts[4], validate=True)
        key = base64.b64decode(parts[5], validate=True)
    except ValueError:
        return None

    is_power_of_two = n > 1 and n & (n - 1) == 0
    if not is_power_of_two or n > MAX_SCRYPT_N or not 1 <= r <= MAX_SCRYPT_R or not 1 <= p <= MAX_SCRYPT_P:
        return None
    if not salt or not key:
        return None
    return n, r, p, salt, key


def _is_legacy(stored_hash: str) -> bool:
    return len(stored_hash) in LEGACY_HEX_LENGTHS and all(c in "0123456789abcdef" for c in stored_hash.lower())


def verify_password(password, stored_hash) -> bool:
    """
    Check a submitted password against the value stored in `guser.pass`.

    Accepts scrypt hashes and the legacy SHA3-256 and MD5 hexdigests. Anything else (including an
    empty or NULL column) never matches.
    """
    if password is None or not stored_hash:
        return False

    parsed = _parse(stored_hash)
    if parsed:
        n, r, p, salt, key = parsed
        return hmac.compare_digest(_scrypt(password, salt, n, r, p, len(key)), key)

    if _is_legacy(stored_hash):
        digest = LEGACY_HEX_LENGTHS[len(stored_hash)](password.encode("utf-8")).hexdigest()
        return hmac.compare_digest(digest, stored_hash.lower())

    return False


def verify_unknown_user(password) -> None:
    """
    Spend the time a real check would take when there is no account to check against, so that the
    response time of a login attempt doesn't show whether an e-mail address is registered.
    """
    _scrypt(password or "", b"\0" * SALT_BYTES, SCRYPT_N, SCRYPT_R, SCRYPT_P, KEY_BYTES)


def needs_rehash(stored_hash) -> bool:
    """True if a stored hash that just verified should be replaced with `hash_password()`'s output."""
    parsed = _parse(stored_hash) if stored_hash else None
    if not parsed:
        return True
    n, r, p, salt, key = parsed
    return (n, r, p) != (SCRYPT_N, SCRYPT_R, SCRYPT_P) or len(salt) < SALT_BYTES or len(key) < KEY_BYTES
