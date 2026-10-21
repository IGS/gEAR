"""
Small lookup and form CGIs: organisms, tags, comments, PubMed citations, workshop event registration,
plus the pure helpers of stats.cgi, get_dataset_comparison.cgi and the legacy process_contact.py.
"""

import importlib.machinery
import importlib.util
import json
import os
import subprocess
import sys
import types
from pathlib import Path
from urllib.parse import urlencode

import pytest
from helpers.cgi_harness import CGI_DIR, run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"


def load_cgi(name, module_name):
    """Import a CGI as a module, restoring sys.stdout if the script redirects it on import."""
    stdout = sys.stdout
    try:
        loader = importlib.machinery.SourceFileLoader(module_name, str(CGI_DIR / name))
        module = importlib.util.module_from_spec(importlib.util.spec_from_loader(loader.name, loader))
        loader.exec_module(module)
    finally:
        if sys.stdout is not stdout:
            sys.stdout.close()
            sys.stdout = stdout
    return module


# ---------------------------------------------------------------- get_organism_list.cgi / get_tag_list.cgi

def test_organism_list_flags_organisms_with_icons():
    db = {"sql": [{"match": "FROM organism", "rows": [[1, "Mouse", "Mus", "musculus", None, 10090],
                                                      [999, "Newt", "Pleurodeles", "waltl", "x", 8319]]}]}
    body = run_cgi("get_organism_list.cgi", db=db).json()
    assert body == {"organisms": [
        {"id": 1, "label": "Mouse", "genus": "Mus", "species": "musculus", "strain": None, "taxon_id": 10090,
         "has_icon": 1},
        {"id": 999, "label": "Newt", "genus": "Pleurodeles", "species": "waltl", "strain": "x", "taxon_id": 8319,
         "has_icon": 0}]}


def test_organism_list_empty():
    assert run_cgi("get_organism_list.cgi").json() == {"organisms": []}


def test_tag_list():
    db = {"sql": [{"match": "FROM tag", "rows": [[1, "bug"], [2, "question"]]}]}
    result = run_cgi("get_tag_list.cgi", db=db)
    assert result.json() == {"success": 1, "tags": [{"id": 1, "label": "bug"}, {"id": 2, "label": "question"}]}
    assert result.writes() == []


# ---------------------------------------------------------------- load_comment.cgi

COMMENT = {"submitter_firstname": "Ann", "submitter_lastname": "Lee", "submitter_email": "ann@example.org",
           "comment_title": "Hi", "comment": "A question"}


def test_comment_without_tags():
    result = run_cgi("load_comment.cgi", form=COMMENT, db={"lastrowid": 42})
    assert result.json() == {"success": 1}
    [insert] = result.writes()
    assert insert["query"].startswith("INSERT INTO comment ")
    assert insert["params"] == ["Ann", "Lee", "ann@example.org", "Hi", "A question", 0]


def test_comment_tags_reuse_existing_and_add_new():
    db = {"lastrowid": 42, "sql": [{"match": "FROM tag", "rows": [["Bug", 7]]}]}
    # Duplicates are dropped and existing tags match case-insensitively
    result = run_cgi("load_comment.cgi", form={**COMMENT, "comment_tag": "bug, new idea, bug"}, db=db)
    assert result.json() == {"success": 1}
    assert [(w["query"].split(" (")[0], w["params"]) for w in result.writes()] == [
        ("INSERT INTO comment", ["Ann", "Lee", "ann@example.org", "Hi", "A question", 0]),
        ("INSERT INTO comment_tag", [7, 42]),
        ("INSERT INTO tag", ["new idea"]),
        ("INSERT INTO comment_tag", [42, 42]),      # the fake's lastrowid is also the new tag's ID
    ]


# ---------------------------------------------------------------- get_citation_from_pubmed_id.cgi

class TestCitation:
    def run(self, tmp_path, pubmed_id="12345", mode="json", body=None):
        log = tmp_path / "requests.jsonl"
        env = {"FAKE_REQUESTS_LOG": str(log), "FAKE_REQUESTS_MODE": mode, "FAKE_REQUESTS_BODY": body or "{}"}
        result = run_cgi("get_citation_from_pubmed_id.cgi", form={"pubmed_id": pubmed_id}, env=env,
                         fake_modules={"requests": FAKES / "cgi_wb_requests.py"})
        calls = [json.loads(line) for line in log.read_text().splitlines()] if log.exists() else []
        return result, calls

    def test_returns_citation(self, tmp_path):
        citation = {"ama": {"orig": "Doe J. A paper."}}
        result, calls = self.run(tmp_path, pubmed_id=" 12345 ", body=json.dumps(citation))
        assert (result.status, result.json()) == (200, citation)
        assert calls == [{"url": "https://pmc.ncbi.nlm.nih.gov/api/ctxp/v1/pubmed/?format=citation&id=12345",
                          "timeout": 10}]

    def test_non_json_answer_is_returned_raw(self, tmp_path):
        result, _ = self.run(tmp_path, mode="text")
        assert result.json() == {"raw": "not json"}

    @pytest.mark.parametrize("mode", ["http_error", "unreachable"])
    def test_upstream_failure_is_bad_gateway(self, tmp_path, mode):
        result, calls = self.run(tmp_path, mode=mode)
        assert result.status == 502 and "error" in result.json()
        assert len(calls) == 1

    def test_missing_id(self, tmp_path):
        result, calls = self.run(tmp_path, pubmed_id="   ")
        assert (result.status, result.json()) == (400, {"error": "No pubmed_id provided"})
        assert calls == []

    @pytest.mark.parametrize("pubmed_id", ["123&format=xml", "123&id=456", "123#x", "../123", "PMID:123", "12 34",
                                           "12a", "-123", "１２３"])
    def test_non_numeric_id_is_refused(self, tmp_path, pubmed_id):
        # Used to be pasted into the upstream URL as-is, so it could add or override query parameters
        result, calls = self.run(tmp_path, pubmed_id=pubmed_id)
        assert (result.status, result.json()) == (400, {"error": "Invalid pubmed_id"})
        assert calls == []


# ---------------------------------------------------------------- workshop event registration

EVENT_ROWS = {"match": "FROM event WHERE", "rows": [[1, "Intro", 2, 1], [2, "Advanced", 10, 0]]}


class TestEventRegistrationList:
    def run(self, registrations, session_id="owner"):
        db = {"sessions": {"owner": 1},
              "sql": [EVENT_ROWS, {"match": "JOIN event_registration", "rows": registrations}]}
        query = {"min_event_id": 1, "max_event_id": 8}
        if session_id:
            query["session_id"] = session_id
        return run_cgi("get_event_registration_list.cgi", db=db, query=query)

    @staticmethod
    def event(attendees, max_attendees, waitlist_size, attending=0, waitlisted=0):
        return {"attendees": attendees, "user_attending": attending, "user_waitlisted": waitlisted,
                "max_attendees": max_attendees, "waitlist_size": waitlist_size}

    def test_anonymous_gets_nothing(self):
        result = self.run([[1, 100, 5]], session_id=None)
        assert result.json() == {}
        assert result.logged("sql") == []

    def test_last_seat_is_attending(self):
        # Event 1 has two seats; the user signed up second (this used to report them as waitlisted)
        result = self.run([[1, 100, 5], [1, 101, 1]])
        assert result.json() == {"1": self.event(2, 2, 1, attending=1), "2": self.event(0, 10, 0)}
        assert "ORDER BY er.id" in result.logged("sql", "JOIN event_registration")[0]["query"]

    def test_past_the_seats_is_waitlisted(self):
        result = self.run([[1, 100, 5], [1, 101, 6], [1, 102, 1], [2, 103, 1]])
        assert result.json() == {"1": self.event(3, 2, 1, waitlisted=1), "2": self.event(1, 10, 0, attending=1)}


class TestSetEventRegistration:
    """Event 1 has 2 seats and a waitlist of 1, as in TestEventRegistrationList."""

    EVENT = {"attendees": 0, "user_attending": 0, "user_waitlisted": 0, "max_attendees": 2, "waitlist_size": 1}

    def run(self, status, session_id="owner", event_id="1", count=1, registered=False, event_exists=True):
        sql = [{"match": "COUNT(er.id)", "rows": [[count]]}]
        if event_exists:
            sql.append({"match": "FROM event WHERE id", "rows": [["Intro", 2, 1]]})
        if registered:
            sql.append({"match": "SELECT id FROM event_registration", "rows": [[5]]})
        form = {"session_id": session_id, "event_id": event_id}
        if status is not None:
            form["registration_status"] = status
        return run_cgi("set_user_event_registration.cgi", db={"sessions": {"owner": 1}, "sql": sql}, form=form)

    @staticmethod
    def assert_refused(result, msg):
        assert result.json() == {"success": 0, "msg": msg}
        assert result.writes() == []
        assert not result.logged("commit")

    def test_register(self):
        result = self.run("1", count=1)
        assert result.json() == {"success": 1, "current_count": 2, "event": self.EVENT,
                                 "msg": "You have successfully registered"}
        assert [(w["query"], w["params"]) for w in result.writes()] == [
            ("INSERT INTO event_registration (event_id, user_id) VALUES (%s, %s)", [1, 1])]
        assert result.logged("commit")
        # The event row is locked while the registrations are counted
        assert result.logged("sql", "FROM event WHERE id = %s FOR UPDATE")

    def test_register_when_seats_are_taken_joins_the_waitlist(self):
        # Like get_event_registration_list.cgi: registrations past max_attendees are waitlisted
        result = self.run("1", count=2)
        body = result.json()
        assert (body["success"], body["current_count"]) == (1, 3)
        assert "waitlist" in body["msg"]
        assert len(result.writes()) == 1 and result.logged("commit")

    def test_register_when_waitlist_is_full_is_refused(self):
        self.assert_refused(self.run("1", count=3), "This event is full")

    def test_register_twice_adds_no_second_registration(self):
        result = self.run("1", count=2, registered=True)
        assert result.json() == {"success": 1, "current_count": 2, "event": self.EVENT,
                                 "msg": "You are already registered"}
        assert result.writes() == []

    def test_unregister(self):
        result = self.run("0")
        assert result.json()["msg"] == "You have successfully unregistered"
        assert [(w["query"], w["params"]) for w in result.writes()] == [
            ("DELETE FROM event_registration WHERE event_id = %s AND user_id = %s", [1, 1])]

    @pytest.mark.parametrize("status", ["7", "", "yes", None])
    def test_unknown_status(self, status):
        # A missing or non-numeric status used to raise before any output (a 500)
        result = self.run(status)
        assert result.returncode == 0, result.stderr
        self.assert_refused(result, "Unrecognized registration status")

    @pytest.mark.parametrize("status", ["1", "0"])
    def test_unknown_event_is_refused(self, status):
        self.assert_refused(self.run(status, event_exists=False), "Unknown event")

    @pytest.mark.parametrize("event_id", ["", "abc", "1; DROP TABLE event"])
    def test_invalid_event_id_is_refused(self, event_id):
        result = self.run("1", event_id=event_id)
        self.assert_refused(result, "Unknown event")
        assert not result.logged("sql", "event_registration")

    @pytest.mark.parametrize("status", ["1", None])
    def test_must_be_logged_in(self, status):
        # Checked before the status is parsed, so a missing status doesn't hide the login message
        result = self.run(status, session_id="nobody")
        assert result.json() == {"success": 0, "msg": "You must be logged in to register or unregister"}
        assert result.writes() == []


# ---------------------------------------------------------------- stats.cgi parsers (no Google Analytics calls)

@pytest.fixture
def stats(monkeypatch):
    """Import stats.cgi, standing in for the Google Analytics client library if it isn't installed."""
    try:
        import google.analytics.data_v1beta  # noqa: F401
    except ImportError:
        data = types.ModuleType("google.analytics.data_v1beta")
        data.BetaAnalyticsDataClient = object
        data.types = types.ModuleType("google.analytics.data_v1beta.types")
        for name in ("DateRange", "Dimension", "Metric", "OrderBy", "RunRealtimeReportRequest", "RunReportRequest"):
            setattr(data.types, name, object)
        try:
            import google
        except ImportError:
            google = types.ModuleType("google")
            monkeypatch.setitem(sys.modules, "google", google)
        analytics = types.ModuleType("google.analytics")
        analytics.data_v1beta = data
        monkeypatch.setattr(google, "analytics", analytics, raising=False)
        monkeypatch.setitem(sys.modules, "google.analytics", analytics)
        monkeypatch.setitem(sys.modules, "google.analytics.data_v1beta", data)
        monkeypatch.setitem(sys.modules, "google.analytics.data_v1beta.types", data.types)
    return load_cgi("stats.cgi", "gear_stats_cgi")


@pytest.mark.parametrize("value, bounds, expected", [
    ("5", {}, 5), ("-3", {}, -3), ("x", {}, 7), (None, {}, 7), ("3.5", {}, 7), ("", {}, 7),
    ("0", {"minimum": 1, "maximum": 365}, 1), ("1000", {"minimum": 1, "maximum": 365}, 365),
    ("42", {"minimum": 1, "maximum": 365}, 42),
])
def test_stats_parse_int(stats, value, bounds, expected):
    assert stats.parse_int(value, 7, **bounds) == expected


@pytest.mark.parametrize("value, expected", [("12", 12), ("12.7", 12), ("1e3", 1000), (3.9, 3), ("abc", 0),
                                             (None, 0), ("", 0)])
def test_stats_safe_int(stats, value, expected):
    assert stats.safe_int(value) == expected


@pytest.mark.parametrize("query, expected", [
    ("", (30, 10, True)),
    ("days=7&top_n=5&realtime=1", (7, 5, True)),
    ("days=0&top_n=100&realtime=No", (1, 50, False)),
    ("days=9999&top_n=0&realtime=false", (365, 1, False)),
    ("days=abc&top_n=&realtime=0", (30, 10, False)),
])
def test_stats_query_params(stats, monkeypatch, query, expected):
    monkeypatch.setenv("QUERY_STRING", query)
    days, top_n, realtime = expected
    assert stats.get_query_params() == {"days": days, "top_n": top_n, "include_realtime": realtime}


# ---------------------------------------------------------------- get_dataset_comparison.cgi helpers

@pytest.fixture(scope="module")
def comparison():
    pytest.importorskip("scanpy")
    return load_cgi("get_dataset_comparison.cgi", "gear_dataset_comparison_cgi")


@pytest.mark.parametrize("x, y, expected", [(4, 2, 2), (2, 4, 2), (3, 3, 1), (3, 0, 3), (0, 3, 3), (0, 0, 0),
                                            (0.5, 2.0, 4.0)])
def test_fold_change_is_larger_over_smaller(comparison, x, y, expected):
    assert comparison.fold_change(x, y) == pytest.approx(expected)


@pytest.mark.parametrize("val, base, expected", [(8, 2, 3), (100, 10, 2), (1, 2, 0), (0, 2, 0), (-1, 2, None)])
def test_get_log(comparison, val, base, expected):
    result = comparison.get_log(val, base)
    assert result == (pytest.approx(expected) if expected is not None else None)


# ---------------------------------------------------------------- process_contact.py (legacy contact form)

def test_process_contact_saves_submission(tmp_path):
    # It writes to ./user-submitted-comments relative to the working directory, so run it from tmp_path
    (tmp_path / "user-submitted-comments").mkdir()
    form = {"submitterName": "Ann", "InputEmail": "ann@example.org", "keep_updated": "on", "comments": "Hello"}
    env = {"PATH": os.environ.get("PATH", ""), "REQUEST_METHOD": "GET", "QUERY_STRING": urlencode(form)}
    proc = subprocess.run([sys.executable, str(CGI_DIR / "process_contact.py")], cwd=tmp_path, env=env,
                          capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.startswith("Content-type: text/html\n\nComment submitted.")
    [saved] = (tmp_path / "user-submitted-comments").iterdir()
    assert saved.read_text() == ("submitter: Ann\ne-mail: ann@example.org\nkeep updated: on\nsec check: None\n"
                                 "comments: Hello\n")
