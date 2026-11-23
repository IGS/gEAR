"""
The unified search API resource, /search (#1383).

The resource is called through Flask's test client. The searchers are replaced with stubs, and
the user and database lookups are monkeypatched, so these tests cover only the request handling:
reading the user and the criteria, running the search, and the response.
"""

import types
from pathlib import Path

import pytest
from flask import Flask
from flask_restful import Api
from gear.search import TYPES, Hit
from resources import search

API_PY = Path(__file__).resolve().parents[2] / "www" / "api" / "api.py"


class StubSearcher:
    """Finds one hit per type and records what each call was given."""

    def __init__(self, kind, calls, error=None):
        self.kind = kind
        self.calls = calls
        self.error = error

    def find(self, criteria, user, cursor):
        self.calls.append(("find", self.kind, criteria, user))
        if self.error:
            raise self.error
        return [Hit(self.kind, f"{self.kind}-1", 1.0, f"A {self.kind}", "2026-01-01")]

    def hydrate(self, keys, user, cursor):
        return {key: {"type": self.kind, "id": key, "title": f"A {self.kind}"} for key in keys}


class FakeConnection:
    def __init__(self):
        self.closed = False
        self.cursor = types.SimpleNamespace(closed=False)

        def close_cursor():
            self.cursor.closed = True
        self.cursor.close = close_cursor

    def get_cursor(self):
        return self.cursor

    def close(self):
        self.closed = True


@pytest.fixture
def api(monkeypatch):
    """A test client for the resource, with .calls (searcher calls) and .connections (database connections)."""
    calls = []
    connections = []
    state = types.SimpleNamespace(error=None)

    def build_searchers():
        return {kind: StubSearcher(kind, calls, state.error) for kind in TYPES}

    def connect():
        connections.append(FakeConnection())
        return connections[-1]

    users = {"owner": types.SimpleNamespace(id=5)}
    monkeypatch.setattr(search, "build_searchers", build_searchers)
    monkeypatch.setattr(search.geardb, "Connection", connect)
    monkeypatch.setattr(search.geardb, "get_user_from_session_id", lambda session_id: users.get(session_id))

    app = Flask(__name__)
    Api(app).add_resource(search.UnifiedSearch, "/search")
    client = app.test_client()
    client.calls = calls
    client.connections = connections
    client.state = state
    return client


def finds(client):
    return [call for call in client.calls if call[0] == "find"]


class TestUnifiedSearch:
    def test_response(self, api):
        resp = api.get("/search?q=hair")
        assert resp.status_code == 200
        body = resp.get_json()
        assert body["success"] == 1 and body["problem"] == ""
        assert body["counts"] == {"all": 6, **{kind: 1 for kind in TYPES}}
        # Every stub hit ties on score and date, so the order is by title; only the types matter here
        assert sorted(result["type"] for result in body["results"]) == sorted(TYPES)
        assert body["pagination"]["total_results"] == 6

    def test_logged_out_searches_as_no_user(self, api):
        api.get("/search")
        assert {user for _, _, _, user in finds(api)} == {None}

    def test_session_cookie_identifies_the_user(self, api):
        api.set_cookie("gear_session_id", "owner")
        api.get("/search")
        assert {user.id for _, _, _, user in finds(api)} == {5}

    def test_unknown_session_is_logged_out(self, api):
        api.set_cookie("gear_session_id", "expired")
        api.get("/search")
        assert {user for _, _, _, user in finds(api)} == {None}

    def test_request_arguments_become_the_criteria(self, api):
        api.get("/search?q=%22hair+cell%22&types=study,ccc&organism_ids=1,2&ownership=yours,group"
                "&date_added=year&sort=title&page=2&limit=5")
        criteria = finds(api)[0][2]
        assert criteria.q == '"hair cell"'
        assert criteria.types == ["study", "ccc"]
        assert criteria.organism_ids == [1, 2]
        assert criteria.ownership == ["yours", "group"]
        assert (criteria.date_added, criteria.sort, criteria.page, criteria.limit) == ("year", "title", 2, 5)

    def test_type_filter_limits_results_but_not_counts(self, api):
        body = api.get("/search?types=study").get_json()
        assert [result["type"] for result in body["results"]] == ["study"]
        assert body["counts"]["all"] == 6

    @pytest.mark.parametrize("query, problem", [
        ("types=nonsense", "Unknown result type: nonsense"),
        ("page=0", "Page must be a whole number greater than 0"),
        ("limit=abc", "Limit must be a whole number greater than 0"),
        ("sort=owner", "Sort must be one of: relevance, title, date"),
        ("organism_ids=1;DROP", "Organism IDs must be numbers"),
    ])
    def test_bad_request_is_a_400_with_the_problem(self, api, query, problem):
        resp = api.get(f"/search?{query}")
        assert resp.status_code == 400
        assert resp.get_json() == {"success": 0, "problem": problem}
        # Rejected before any database work
        assert api.connections == [] and api.calls == []

    def test_failed_search_is_a_500_without_details(self, api):
        api.state.error = RuntimeError("SQL syntax error near 'secret'")
        resp = api.get("/search?q=hair")
        assert resp.status_code == 500
        assert resp.get_json() == {"success": 0,
                                    "problem": "The search could not be run. Please try different search terms."}

    @pytest.mark.parametrize("error", [None, RuntimeError("boom")])
    def test_database_connection_is_closed(self, api, error):
        api.state.error = error
        api.get("/search")
        [connection] = api.connections
        assert connection.closed and connection.cursor.closed


class TestBuildSearchers:
    def test_one_searcher_per_type(self, monkeypatch):
        monkeypatch.setattr(search.geardb, "domain_label", "nemo")
        searchers = search.build_searchers()
        assert set(searchers) == set(TYPES)
        for searcher in searchers.values():
            assert callable(searcher.find) and callable(searcher.hydrate)
        # Study pages are filtered to this site's portal
        assert searchers["study"].domain_label == "nemo"

    def test_datasets_find_nothing_until_their_searcher_is_written(self):
        # TODO #1383: replace when the dataset searcher lands (after the domain-specific landing page merge)
        assert search.build_searchers()["dataset"].find(search.parse_criteria({}), None, None) == []


def test_route_is_registered():
    # api.py can't be imported in tests (it sets memory limits and queries the database at import)
    source = API_PY.read_text()
    assert "from resources.search import UnifiedSearch" in source
    assert "api.add_resource(UnifiedSearch, '/search')" in source
