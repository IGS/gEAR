"""
fake_geardb.py - A stand-in for lib/geardb.py used when running CGI scripts in tests.

The fake is configured by a JSON "spec" file whose path is in the GEAR_FAKE_DB environment
variable, and it appends a JSON line for every SQL statement, commit and object save to the
file named by GEAR_FAKE_DB_LOG, so tests can assert what was (or was not) written.

Spec keys (all optional):
    sessions        {session_id: user_id}
    users           {user_id: {email, user_name, ...}}
    datasets        [{id, share_id, owner_id, is_downloadable, is_public, dtype, title, ...,
                      archive_path, h5ad_path}]
    layouts         [{id, share_id, user_id, label, is_public}]
    gene_carts      [{id, share_id, user_id, gctype, label}]
    metadata        {share_id: csv text}
    verification_short   value returned by get_verification_code_short_form()
    sql             [{"match": "substring", "rows": [[...], ...] | "echo"}]
                    First matching rule wins. "echo" returns [[first param]].
    fail_sql        ["substring", ...]  -- execute() raises for matching statements
    lastrowid       value of cursor.lastrowid (default 1)

Queries against user_session are answered from "sessions" unless an sql rule matches first.
Any geardb attribute the fake doesn't provide raises AttributeError naming this file, so a
missing stub is obvious rather than silently wrong.
"""

import json
import os
import types
from pathlib import Path

# fail_sql raises a real mysql.connector.Error, so CGIs that catch only database errors handle it
try:
    from mysql.connector import Error as _DatabaseError
except ImportError:
    _DatabaseError = RuntimeError

_SPEC = json.loads(Path(os.environ["GEAR_FAKE_DB"]).read_text()) if os.environ.get("GEAR_FAKE_DB") else {}
_LOG = os.environ.get("GEAR_FAKE_DB_LOG")

# Like the real geardb, make json.dumps() accept gEAR objects: use a class's _serialize_json()
#  if it has one, else fall back to str() (production CGIs rely on this)
def _json_default(self, obj):
    try:
        return getattr(obj.__class__, "_serialize_json", _json_default.default)(obj)
    except Exception:
        return str(obj)


_json_default.default = json.JSONEncoder().default
json.JSONEncoder.default = _json_default

servercfg = {}
analysis_base_dir = "/tmp"
domain_url = "http://localhost"
domain_label = "gEAR"


def _log(event, **fields):
    """Record an event (SQL, commit, save) for the test to inspect."""
    if _LOG:
        with open(_LOG, "a") as fh:
            fh.write(json.dumps({"event": event, **fields}, default=str) + "\n")


class _User:
    def __init__(self, user_id):
        self.id = int(user_id)
        info = _SPEC.get("users", {}).get(str(user_id), {})
        self.email = info.get("email", f"user{user_id}@example.org")
        self.user_name = info.get("user_name", f"User {user_id}")
        self.is_admin = info.get("is_admin", 0)
        self.help_id = info.get("help_id", f"help-{user_id}")


def get_user_from_session_id(session_id=None):
    user_id = _SPEC.get("sessions", {}).get(session_id or "")
    return _User(user_id) if user_id is not None else None


def get_user_id_from_session_id(session_id=None):
    user = get_user_from_session_id(session_id)
    return user.id if user else None


# ---------------------------------------------------------------- layouts (collections)

class _Layout:
    def __init__(self, spec):
        self.id = spec.get("id", 1)
        self.share_id = spec.get("share_id")
        self.user_id = spec.get("user_id")
        self.label = spec.get("label", "A collection")
        self.is_public = spec.get("is_public", 0)
        self.members = []

    def load(self):
        pass

    def get_members(self):
        return self.members

    def remove_member_by_dataset_id(self, dataset_id):
        _log("layout.remove_member_by_dataset_id", layout=self.share_id, dataset_id=dataset_id)

    def remove_member_by_display_id(self, display_id):
        _log("layout.remove_member_by_display_id", layout=self.share_id, display_id=display_id)

    def remove(self):
        _log("layout.remove", layout=self.share_id)

    def save(self):
        _log("layout.save", layout=self.share_id, label=self.label)

    def save_change(self, attribute, value):
        _log("layout.save_change", layout=self.share_id, attribute=attribute, value=value)


def get_layout_by_share_id(layout_share_id):
    for spec in _SPEC.get("layouts", []):
        if spec.get("share_id") == layout_share_id:
            return _Layout(spec)
    return None


# ---------------------------------------------------------------- datasets

class Dataset(types.SimpleNamespace):
    """Fake geardb.Dataset built from a spec entry (or constructed directly by a CGI)."""

    def __init__(self, **kwargs):
        spec = next((d for d in _SPEC.get("datasets", []) if d.get("id") == kwargs.get("id")), {})
        defaults = {"is_downloadable": 1, "is_public": 1, "dtype": "single-cell-rnaseq", "title": "A dataset",
                    "ldesc": None, "pubmed_id": None, "geo_id": None, "organism_id": 1, "owner_id": None,
                    "share_id": None, "schematic_image": None, "load_status": "completed"}
        merged = {**defaults, **spec, **kwargs}
        super().__init__(**merged)

    def get_source_archive_path(self):
        return getattr(self, "archive_path", None)

    def get_tarball_path(self):
        return getattr(self, "archive_path", None) or f"/nonexistent/{self.id}.tar.gz"

    def get_file_path(self, *args, **kwargs):
        return getattr(self, "h5ad_path", None) or f"/nonexistent/{self.id}.h5ad"

    def save_change(self, attribute, value):
        _log("dataset.save_change", dataset=self.id, attribute=attribute, value=value)

    def _serialize_json(self):
        return dict(vars(self))


def get_dataset_by_id(d_id=None, include_shape=None, **kwargs):
    d_id = d_id or kwargs.get("dataset_id")
    for spec in _SPEC.get("datasets", []):
        if spec.get("id") == d_id:
            return Dataset(**spec)
    return None


def get_dataset_by_share_id(share_id, *args, **kwargs):
    for spec in _SPEC.get("datasets", []):
        if spec.get("share_id") == share_id:
            return Dataset(**spec)
    return None


class DatasetCollection:
    """Returns no datasets; tests of search CGIs assert on the logged SQL instead."""

    def __init__(self, *args, **kwargs):
        self.datasets = []

    def get_by_dataset_ids(self, ids=None, **kwargs):
        return self.datasets


def get_metadata_by_share_id(share_id):
    return _SPEC.get("metadata", {}).get(share_id)


def get_displays_by_user_id(user_id=None, dataset_id=None, **kwargs):
    return []


def get_organism_id_by_taxon_id(taxon_id):
    return _SPEC.get("organism_id", 1)


def get_verification_code_short_form(long_form):
    return _SPEC.get("verification_short", "SHORT")


# ---------------------------------------------------------------- gene carts (lists)

def get_gene_cart_by_share_id(share_id, *args, **kwargs):
    for spec in _SPEC.get("gene_carts", []):
        if spec.get("share_id") == share_id:
            return types.SimpleNamespace(**spec)
    return None


# ---------------------------------------------------------------- database connection

class _Cursor:
    def __init__(self):
        self.rows = []
        self.lastrowid = _SPEC.get("lastrowid", 1)
        self.rowcount = 0

    def execute(self, query, params=()):
        flat = " ".join(str(query).split())
        params = list(params or ())
        _log("sql", query=flat, params=params)
        for pattern in _SPEC.get("fail_sql", []):
            if pattern in flat:
                raise _DatabaseError(f"fake database error for: {pattern}")
        self.rows = []
        for rule in _SPEC.get("sql", []):
            if rule["match"] in flat:
                self.rows = [[params[0]]] if rule["rows"] == "echo" else [list(r) for r in rule["rows"]]
                break
        else:
            if "FROM user_session" in flat and params:
                user_id = _SPEC.get("sessions", {}).get(params[0])
                self.rows = [[user_id]] if user_id is not None else []
        self.rowcount = len(self.rows)

    def __iter__(self):
        return iter([tuple(r) for r in self.rows])

    def fetchone(self):
        return tuple(self.rows[0]) if self.rows else None

    def fetchall(self):
        return [tuple(r) for r in self.rows]

    def close(self):
        pass


class Connection:
    def __init__(self):
        pass

    def get_cursor(self, use_dict=False, **kwargs):
        return _Cursor()

    def commit(self):
        _log("commit")

    def close(self):
        pass


def __getattr__(name):
    raise AttributeError(
        f"fake geardb has no attribute {name!r}; add a stub for it in tests/python/helpers/fake_geardb.py"
    )
