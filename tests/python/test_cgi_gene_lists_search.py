"""Gene list and search CGIs: error handling, ownership checks and exactly one JSON response."""

from pathlib import Path

import pytest

from helpers.cgi_harness import run_cgi

FAKES = Path(__file__).resolve().parent / "fakes"
USERHISTORY = {"gear.userhistory": FAKES / "accounts_userhistory.py"}

# Session "owner" is user 1, who owns gene list GC1 (id 5); "other" is user 2.
SESSIONS = {"owner": 1, "other": 2}
UNWEIGHTED = {"id": 5, "share_id": "GC1", "user_id": 1, "gctype": "unweighted-list", "label": "My list"}


class TestGetGeneCartMembers:
    def test_missing_share_id(self):
        body = run_cgi("get_gene_cart_members.cgi").json()
        assert body["gene_symbols"] == [] and body["success"] == 0
        assert "missing gene cart share ID" in body["error"]

    def test_unknown_cart(self):
        body = run_cgi("get_gene_cart_members.cgi", query={"share_id": "NOPE"}).json()
        assert body["success"] == 0 and "NOPE" in body["error"]

    def test_unweighted_cart(self):
        db = {"gene_carts": [UNWEIGHTED],
              "sql": [{"match": "FROM gene_cart_member", "rows": [[11, "Sox2"], [12, "Pax6"]]}]}
        result = run_cgi("get_gene_cart_members.cgi", query={"share_id": "GC1"}, db=db)
        assert result.json() == {
            "gene_symbols": [{"id": 11, "label": "Sox2"}, {"id": 12, "label": "Pax6"}],
            "success": 1,
        }
        (query,) = result.logged("sql", "FROM gene_cart_member")
        assert query["params"] == [5]

    def test_weighted_cart(self, tmp_cart):
        share_id = tmp_cart("id\tgene\tP1\nE1\tSox2\t0.5\n")
        db = {"gene_carts": [{"id": 6, "share_id": share_id, "user_id": 1, "gctype": "weighted-list"}]}
        body = run_cgi("get_gene_cart_members.cgi", query={"share_id": share_id}, db=db).json()
        assert body == {"gene_symbols": [{"id": "E1", "label": "Sox2"}], "success": 1}

    def test_unknown_type_is_one_error(self):
        db = {"gene_carts": [{**UNWEIGHTED, "gctype": "mystery"}]}
        body = run_cgi("get_gene_cart_members.cgi", query={"share_id": "GC1"}, db=db).json()
        assert body["success"] == 0 and "mystery" in body["error"]


class TestRemoveGeneCart:
    # Ownership is checked with its own query: SELECT id, user_id FROM gene_cart WHERE id = %s
    DB = {"sessions": SESSIONS, "gene_carts": [UNWEIGHTED],
          "sql": [{"match": "FROM gene_cart WHERE", "rows": [[5, 1]]}]}

    @pytest.mark.parametrize("session", ["other", None])
    def test_non_owner_is_refused(self, session):
        query = {"share_id": "GC1"}
        if session:
            query["session_id"] = session
        result = run_cgi("remove_gene_cart.cgi", db=self.DB, query=query)
        body = result.json()
        assert body["success"] == 0 and "gene list" in body["error"]
        assert result.writes() == []

    def test_unknown_cart(self):
        body = run_cgi("remove_gene_cart.cgi", db=self.DB, query={"session_id": "owner", "share_id": "NOPE"}).json()
        assert body == {"success": 0, "error": "Invalid share_id."}


@pytest.mark.parametrize("script", ["search_gene_carts.cgi", "search_datasets.cgi"])
@pytest.mark.parametrize("query, problem", [
    ({"page": "abc"}, "Page must be a number"),
    ({"page": "0"}, "Page must be greater than 0"),
    ({"limit": "x"}, "Limit must be a number"),
    ({"limit": "0"}, "Limit must be greater than 0"),
])
def test_search_rejects_bad_paging(script, query, problem):
    result = run_cgi(script, query=query, fake_modules=USERHISTORY)
    body = result.json()     # exactly one JSON document
    assert body["success"] == 0
    assert body["problem"] == problem
    assert not result.logged("sql")


def test_get_metadata_from_geo_without_id():
    result = run_cgi("get_metadata_from_geo.cgi")
    assert result.json() == {}


def test_get_metadata_from_geo_unrecognised_id():
    # Neither GSE nor GSM: no lookup is attempted
    result = run_cgi("get_metadata_from_geo.cgi", query={"geo_id": "XYZ123"})
    assert result.json() == {}
