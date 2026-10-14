"""
cgi_harness.py - Run gEAR CGI scripts in a subprocess against a fake database.

Example:
    result = run_cgi("rename_layout.cgi",
                     query={"session_id": "owner", "layout_share_id": "L1", "layout_name": "New"},
                     db={"sessions": {"owner": 1}, "layouts": [{"share_id": "L1", "user_id": 1}]})
    assert result.json()["layout_label"] == "New"
    assert result.logged("layout.save")

Each run gets its own process, so scripts that swap sys.stdout or call sys.exit() can't affect
the test process. The CGI runs from www/cgi (many scripts find lib/ relative to that directory).
"""

import json
import os
import subprocess
import sys
import tempfile
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urlencode

HELPERS_DIR = Path(__file__).resolve().parent
REPO_ROOT = HELPERS_DIR.parents[2]
CGI_DIR = REPO_ROOT / "www" / "cgi"


@dataclass
class CGIResult:
    """Parsed output of one CGI run."""

    returncode: int
    headers: dict
    body: bytes
    stderr: str
    db_log: list = field(default_factory=list)

    @property
    def status(self) -> int:
        """HTTP status code (200 unless the script printed a Status header)."""
        return int(self.headers.get("status", "200").split()[0])

    @property
    def text(self) -> str:
        return self.body.decode("utf-8", errors="replace")

    def json_documents(self) -> list:
        """All JSON documents in the body, in order (a well-behaved CGI prints exactly one)."""
        decoder = json.JSONDecoder()
        text, pos, docs = self.text, 0, []
        while True:
            while pos < len(text) and text[pos].isspace():
                pos += 1
            if pos >= len(text):
                return docs
            doc, pos = decoder.raw_decode(text, pos)
            docs.append(doc)

    def json(self):
        """The single JSON document in the body (fails if there are zero or several)."""
        docs = self.json_documents()
        assert len(docs) == 1, f"expected one JSON document, got {len(docs)}: {self.text[:500]!r}"
        return docs[0]

    def logged(self, event: str, contains: str = "") -> list:
        """Fake-database log entries of the given event whose query/fields contain `contains`."""
        return [e for e in self.db_log if e["event"] == event and contains in json.dumps(e)]

    def writes(self) -> list:
        """SQL statements that modify data (INSERT/UPDATE/DELETE)."""
        return [e for e in self.logged("sql") if e["query"].split()[0].upper() in ("INSERT", "UPDATE", "DELETE")]


def _parse_output(raw: bytes):
    """Split CGI output into headers (lower-cased names) and body at the first blank line."""
    normalized = raw.replace(b"\r\n", b"\n")
    head, sep, body = normalized.partition(b"\n\n")
    if not sep:
        return {}, normalized
    headers = {}
    for line in head.decode("latin-1").splitlines():
        if ":" in line:
            name, value = line.split(":", 1)
            headers[name.strip().lower()] = value.strip()
    return headers, body


def _multipart(form: dict, files: dict):
    """Build a multipart/form-data body. files maps field name -> (filename, bytes)."""
    boundary = uuid.uuid4().hex
    parts = []
    for name, value in (form or {}).items():
        parts.append(f'--{boundary}\r\nContent-Disposition: form-data; name="{name}"\r\n\r\n{value}\r\n'.encode())
    for name, (filename, content) in (files or {}).items():
        parts.append(
            f'--{boundary}\r\nContent-Disposition: form-data; name="{name}"; filename="{filename}"\r\n'
            f"Content-Type: application/octet-stream\r\n\r\n".encode() + content + b"\r\n"
        )
    parts.append(f"--{boundary}--\r\n".encode())
    return b"".join(parts), f"multipart/form-data; boundary={boundary}"


def run_cgi(script, query=None, form=None, files=None, cookies=None, db=None, fake_modules=None,
            env=None, timeout=120) -> CGIResult:
    """
    Run www/cgi/<script> and return its parsed output.

    query: dict sent as QUERY_STRING (GET). form/files: sent as a multipart POST body.
    cookies: dict sent as HTTP_COOKIE. db: fake database spec (see fake_geardb.py).
    fake_modules: {module name: path to .py} to replace, e.g. {"gear.analysis": ".../fake_analysis.py"}.
    """
    with tempfile.TemporaryDirectory() as tmp:
        spec_path = Path(tmp) / "db.json"
        log_path = Path(tmp) / "db_log.jsonl"
        spec_path.write_text(json.dumps(db or {}))

        run_env = {
            "PATH": os.environ.get("PATH", ""),
            "HOME": os.environ.get("HOME", tmp),
            "GEAR_FAKE_DB": str(spec_path),
            "GEAR_FAKE_DB_LOG": str(log_path),
            "GEAR_FAKE_MODULES": json.dumps({k: str(v) for k, v in (fake_modules or {}).items()}),
            "MPLBACKEND": "Agg",
            "NUMBA_CACHE_DIR": tmp,
        }
        stdin = b""
        if form is not None or files is not None:
            stdin, content_type = _multipart(form, files)
            run_env.update(REQUEST_METHOD="POST", CONTENT_TYPE=content_type, CONTENT_LENGTH=str(len(stdin)))
        else:
            run_env.update(REQUEST_METHOD="GET", QUERY_STRING=urlencode(query or {}))
        if cookies:
            run_env["HTTP_COOKIE"] = "; ".join(f"{k}={v}" for k, v in cookies.items())
        run_env.update(env or {})

        proc = subprocess.run(
            [sys.executable, str(HELPERS_DIR / "run_cgi.py"), str(CGI_DIR / script)],
            input=stdin, capture_output=True, cwd=CGI_DIR, env=run_env, timeout=timeout,
        )
        headers, body = _parse_output(proc.stdout)
        db_log = [json.loads(line) for line in log_path.read_text().splitlines()] if log_path.exists() else []
        return CGIResult(proc.returncode, headers, body, proc.stderr.decode(errors="replace"), db_log)
