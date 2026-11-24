"""Gene list and search CGIs: error handling, ownership checks and exactly one JSON response."""

from pathlib import Path
from types import SimpleNamespace

import pytest
from gear.search import access
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


SEARCH_SCRIPTS = pytest.mark.parametrize("script", ["search_gene_carts.cgi", "search_datasets.cgi"])


def match_queries(result):
    return [e for e in result.logged("sql") if "MATCH(" in e["query"]]


@SEARCH_SCRIPTS
def test_search_text_is_converted_for_boolean_mode(script):
    result = run_cgi(script, query={"search_terms": "Garcia-Añoveros)", "sort_by": "relevance"},
                     fake_modules=USERHISTORY)
    body = result.json()
    assert body["success"] == 1, body
    # search_gene_carts.cgi only runs its count query when the main query found something
    main_query, *count_queries = match_queries(result)
    for entry in (main_query, *count_queries):
        assert "AGAINST(%s IN BOOLEAN MODE)" in entry["query"]
        assert '"%s"' not in entry["query"]
        assert "Garcia-Añoveros)" not in str(entry["params"])
    # SELECT (relevance score) and WHERE both get the converted text
    assert main_query["params"].count('"Garcia Añoveros"') == 2
    for entry in count_queries:
        assert entry["params"].count('"Garcia Añoveros"') == 1


@pytest.mark.parametrize("script, like_column", [("search_gene_carts.cgi", "gc.label"),
                                                  ("search_datasets.cgi", "d.title")])
@pytest.mark.parametrize("search_terms, fulltext, likes", [
    # #636: "Tono" finds "Tonotopic" (prefix), and words inside underscore names are found with LIKE.
    #  Several unquoted words also rank rows containing them as a phrase first (>"...").
    ("Tono Groves_potential ab", 'Tono* Groves_potential* ab >"Tono Groves_potential ab"',
     ["%Tono%", "%Groves\\_potential%"]),
    # Quoted text matches only as that phrase, so "cell" alone doesn't match every dataset
    ('"hair cell" cochlea', '"hair cell" cochlea*', ["%hair_cell%", "%cochlea%"]),
])
def test_search_matches_partial_words_and_phrases(script, like_column, search_terms, fulltext, likes):
    result = run_cgi(script, query={"search_terms": search_terms, "sort_by": "relevance"},
                     fake_modules=USERHISTORY)
    assert result.json()["success"] == 1
    main_query, *count_queries = match_queries(result)
    assert f"OR {like_column} LIKE %s OR {like_column} LIKE %s)" in main_query["query"]
    # SELECT (relevance score) gets only the full-text query; WHERE gets it followed by the LIKE patterns
    params = main_query["params"]
    assert params[0] == fulltext and params.count(fulltext) == 2
    where = params.index(fulltext, 1)
    assert params[where + 1:where + 3] == likes
    for entry in count_queries:
        where = entry["params"].index(fulltext)
        assert entry["params"][where + 1:where + 3] == likes


@SEARCH_SCRIPTS
def test_punctuation_only_search_lists_without_matching(script):
    result = run_cgi(script, query={"search_terms": "()", "sort_by": "relevance"}, fake_modules=USERHISTORY)
    assert result.json()["success"] == 1
    assert match_queries(result) == []
    assert result.logged("sql")


@SEARCH_SCRIPTS
def test_failed_query_returns_one_json_problem(script):
    result = run_cgi(script, query={"search_terms": "hair cell"}, fake_modules=USERHISTORY,
                     db={"fail_sql": ["MATCH("]})
    assert result.returncode == 0, result.stderr[-2000:]
    body = result.json()     # exactly one JSON document
    assert body["success"] == 0
    assert body["problem"].startswith("The search could not be run")


class TestGetDatasetListScopeSearch:
    """get_dataset_list.cgi's legacy "scope" search uses the same boolean-mode conversion."""

    def run(self, search_terms):
        return run_cgi("get_dataset_list.cgi", db={"sessions": SESSIONS},
                       query={"session_id": "owner", "scope": "others", "search_terms": search_terms})

    def test_punctuated_terms_match_as_a_phrase(self):
        result = self.run("Garcia-Añoveros)")
        assert result.returncode == 0, result.stderr[-2000:]
        assert result.json() == {"datasets": []}
        (search,) = match_queries(result)
        assert search["params"] == [1, '"Garcia Añoveros"', '"Garcia Añoveros"']
        # Tags are still looked up by the terms as typed
        assert result.logged("sql", "dataset_tag")[0]["params"] == ["Garcia-Añoveros)"]

    def test_punctuation_only_lists_without_matching(self):
        result = self.run("()")
        assert result.returncode == 0, result.stderr[-2000:]
        assert match_queries(result) == []
        assert not result.logged("sql", "dataset_tag")


@pytest.mark.parametrize("script, table, alias", [("search_datasets.cgi", "dataset d", "d"),
                                                  ("search_gene_carts.cgi", "gene_cart gc", "gc")])
class TestSearchByShareId:
    """
    #1383: the unified search links to one dataset in Dataset Explorer, or one gene list in Gene List
    Manager, by share_id. The search then returns only that item, if the user may see it.
    """

    def main_query(self, result, table):
        (entry,) = [e for e in result.logged("sql") if f"FROM {table}" in e["query"] and "COUNT(" not in e["query"]]
        return entry

    def test_share_id_narrows_to_that_item(self, script, table, alias):
        result = run_cgi(script, fake_modules=USERHISTORY, query={"share_id": "ABC123"})
        assert result.json()["success"] == 1
        entry = self.main_query(result, table)
        assert f"{alias}.share_id = %s" in entry["query"]
        assert "ABC123" in entry["params"]

    def test_visibility_still_applies(self, script, table, alias):
        # Logged out: still public only
        result = run_cgi(script, fake_modules=USERHISTORY, query={"share_id": "ABC123"})
        assert f"{alias}.is_public = 1" in self.main_query(result, table)["query"]

    def test_other_filters_are_ignored(self, script, table, alias):
        # The page's remembered filters must not hide the linked item, so the search uses only the
        #  share_id and the default visibility (which includes public items, unlike ownership=yours)
        result = run_cgi(script, db={"sessions": SESSIONS}, fake_modules=USERHISTORY,
                         query={"session_id": "owner", "share_id": "ABC123", "search_terms": "hair",
                                "ownership": "yours", "organism_ids": "1", "date_added": "week"})
        assert result.json()["success"] == 1
        entry = self.main_query(result, table)
        query = entry["query"]
        assert "MATCH(" not in query and "LIKE" not in query
        assert "organism_id in" not in query.lower() and "INTERVAL" not in query
        assert f"{alias}.is_public = 1" in query
        assert "hair" not in str(entry["params"])


@pytest.mark.parametrize("script, table, visibility", [
    ("search_datasets.cgi", "dataset d", access.dataset_visibility),
    ("search_gene_carts.cgi", "gene_cart gc", access.gene_cart_visibility),
])
@pytest.mark.parametrize("session_id, ownership", [
    (None, ""), (None, "yours"), ("owner", ""), ("owner", "yours,group"), ("owner", "shared"),
    ("owner", "public"), ("owner", "nonsense"),
])
def test_search_uses_the_shared_visibility_rules(script, table, visibility, session_id, ownership):
    """#1383: these pages and the unified search take who-sees-what from one place, gear.search.access."""
    query = {"ownership": ownership, **({"session_id": session_id} if session_id else {})}
    result = run_cgi(script, db={"sessions": SESSIONS}, fake_modules=USERHISTORY, query=query)
    assert result.json()["success"] == 1
    (entry,) = [e for e in result.logged("sql") if f"FROM {table}" in e["query"] and "COUNT(" not in e["query"]]

    user = SimpleNamespace(id=SESSIONS[session_id]) if session_id else None
    clause, params = visibility(user, ownership.split(",") if ownership else [])
    assert clause in entry["query"]
    assert entry["params"] == params


def test_dataset_search_by_share_id_counts_only_that_dataset():
    """#1383: the result count (and so the pagination) uses the same share_id condition."""
    db = {"sql": [{"match": "SELECT d.id", "rows": [["DS1"]]}, {"match": "SELECT COUNT(*)", "rows": [[1]]}]}
    result = run_cgi("search_datasets.cgi", db=db, fake_modules=USERHISTORY, query={"share_id": "ABC123"})
    assert result.json()["pagination"]["total_results"] == 1
    (count,) = result.logged("sql", "COUNT(*)")
    assert "d.share_id = %s" in count["query"] and "ABC123" in count["params"]
