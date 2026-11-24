"""
Unified search library (#1383): visibility clauses, request parsing, and how results from the
per-type searchers are counted, filtered, sorted, paged and hydrated.

Searchers are replaced by in-memory stubs, so no database is needed.
"""

import types

import pytest
from gear import search
from gear.search import Criteria, Hit, access, parse_criteria, run_search

OWNER = types.SimpleNamespace(id=7)


def placeholders_match(clause, params):
    return clause.count("%s") == len(params)


# ---------------------------------------------------------------- visibility clauses


class TestDatasetVisibility:
    def test_logged_out_sees_only_public_even_with_ownership_filter(self):
        clause, params = access.dataset_visibility(None, ["yours", "shared"])
        assert clause == "(d.is_public = 1)"
        assert params == []

    def test_logged_in_default_is_public_owned_or_shared(self):
        clause, params = access.dataset_visibility(OWNER, [])
        assert "d.is_public = 1" in clause
        assert "d.owner_id = %s" in clause
        assert "dataset_shares" in clause and "is_allowed = 1" in clause
        assert "user_group_membership" not in clause
        assert params == [7, 7]
        assert placeholders_match(clause, params)

    @pytest.mark.parametrize("ownership, must_contain, must_not_contain", [
        (["yours"], "d.owner_id = %s", "d.is_public"),
        (["public"], "d.is_public = 1", "d.owner_id"),
        (["shared"], "dataset_shares", "d.is_public"),
        (["group"], "user_group_membership", "d.is_public"),
    ])
    def test_each_ownership_filter_alone(self, ownership, must_contain, must_not_contain):
        clause, params = access.dataset_visibility(OWNER, ownership)
        assert must_contain in clause
        assert must_not_contain not in clause
        assert placeholders_match(clause, params)

    def test_ownership_filters_combine_with_or(self):
        clause, params = access.dataset_visibility(OWNER, ["yours", "public"])
        assert clause == "(d.owner_id = %s OR d.is_public = 1)"
        assert params == [7]

    def test_unrecognized_ownership_matches_nothing(self):
        assert access.dataset_visibility(OWNER, ["bogus"]) == ("(0 = 1)", [])


class TestGeneCartVisibility:
    def test_logged_out_sees_only_public(self):
        assert access.gene_cart_visibility(None, ["yours"]) == ("(gc.is_public = 1)", [])

    def test_logged_in_default_is_public_owned_or_group(self):
        clause, params = access.gene_cart_visibility(OWNER, [])
        assert "gc.is_public = 1" in clause
        assert "gc.user_id = %s" in clause
        assert "user_group_membership" in clause
        assert placeholders_match(clause, params)
        assert set(params) == {7}

    def test_shared_alone_matches_nothing(self):
        # Sharing gene lists with a user isn't implemented (same as search_gene_carts.cgi)
        assert access.gene_cart_visibility(OWNER, ["shared"]) == ("(0 = 1)", [])


class TestLayoutVisibility:
    def test_logged_out_sees_public_and_domain(self):
        assert access.layout_visibility(None, []) == ("(l.is_public = 1 OR l.is_domain = 1)", [])

    def test_logged_in_default_adds_owned_and_group(self):
        clause, params = access.layout_visibility(OWNER, [])
        assert "l.is_public = 1" in clause and "l.is_domain = 1" in clause
        assert "l.user_id = %s" in clause
        assert "layout_group_membership" in clause
        assert placeholders_match(clause, params)
        assert set(params) == {7}

    def test_public_filter_includes_domain_collections(self):
        assert access.layout_visibility(OWNER, ["public"]) == ("(l.is_public = 1 OR l.is_domain = 1)", [])

    def test_shared_alone_matches_nothing(self):
        # Collections are shared by link only, so there is no "shared with you" list to search
        assert access.layout_visibility(OWNER, ["shared"]) == ("(0 = 1)", [])


# ---------------------------------------------------------------- request parsing


class TestParseCriteria:
    def test_defaults(self):
        assert parse_criteria({}) == Criteria()

    def test_values_are_parsed(self):
        criteria = parse_criteria({
            "q": "  hair cell ", "types": "genelist,ccc", "organism_ids": "1,3",
            "ownership": "yours,public", "date_added": "month", "sort": "title", "page": "2", "limit": "10",
        })
        assert criteria == Criteria(q="hair cell", types=["genelist", "ccc"], organism_ids=[1, 3],
                                    ownership=["yours", "public"], date_added="month", sort="title",
                                    page=2, limit=10)

    def test_selecting_every_type_means_all(self):
        assert parse_criteria({"types": ",".join(search.TYPES)}).types == []

    def test_limit_is_capped(self):
        assert parse_criteria({"limit": "5000"}).limit == search.MAX_LIMIT

    @pytest.mark.parametrize("args", [
        {"types": "dataset,bogus"},
        {"page": "0"},
        {"page": "x"},
        {"limit": "0"},
        {"organism_ids": "1,abc"},
        {"date_added": "decade"},
        {"sort": "owner"},
    ])
    def test_bad_values_raise(self, args):
        with pytest.raises(ValueError):
            parse_criteria(args)


# ---------------------------------------------------------------- combining searchers


class StubSearcher:
    """Returns fixed hits; hydrate() records the keys it was asked for."""

    def __init__(self, hits=(), missing=()):
        self.hits = list(hits)
        self.missing = set(missing)     # keys that vanished between find() and hydrate()
        self.hydrated = []

    def find(self, criteria, user, cursor):
        return list(self.hits)

    def hydrate(self, keys, user, cursor):
        self.hydrated.append(list(keys))
        return {key: {"key": key} for key in keys if key not in self.missing}


def stub_searchers(**hits_by_type):
    return {kind: StubSearcher(hits_by_type.get(kind, ())) for kind in search.TYPES}


def keys(response):
    return [result["key"] for result in response["results"]]


class TestRunSearch:
    def test_counts_cover_every_type_and_ignore_the_type_filter(self):
        searchers = stub_searchers(
            dataset=[Hit("dataset", "d1"), Hit("dataset", "d2")],
            genelist=[Hit("genelist", "g1")],
            ccc=[Hit("ccc", "c1")],
        )
        response = run_search(Criteria(types=["genelist"]), None, None, searchers)
        assert response["counts"] == {"all": 4, "dataset": 2, "collection": 0, "genelist": 1,
                                      "study": 0, "projection": 0, "ccc": 1}
        assert keys(response) == ["g1"]

    def test_type_filters_combine(self):
        searchers = stub_searchers(
            dataset=[Hit("dataset", "d1", date="2026-01-01")],
            genelist=[Hit("genelist", "g1", date="2026-02-01")],
            ccc=[Hit("ccc", "c1", date="2026-03-01")],
        )
        response = run_search(Criteria(types=["dataset", "ccc"], sort="date"), None, None, searchers)
        assert keys(response) == ["c1", "d1"]

    def test_relevance_normalises_scores_within_each_type(self):
        # Raw full-text scores from different tables aren't comparable, so each is scaled by
        #  the best score of its own type. Equal scaled scores go newest first.
        searchers = stub_searchers(
            dataset=[Hit("dataset", "A", 10, date="2026-01-01"), Hit("dataset", "B", 5, date="2026-01-01")],
            genelist=[Hit("genelist", "X", 0.5, date="2026-05-01"), Hit("genelist", "Y", 0.1, date="2026-01-01")],
        )
        response = run_search(Criteria(q="hair"), None, None, searchers)
        assert keys(response) == ["X", "A", "B", "Y"]

    def test_relevance_ties_take_turns_between_types(self):
        # Many hits of one type can share the top scaled score. Ties take turns between types
        #  (newest first within each turn), so one recent type can't fill the first page.
        searchers = stub_searchers(
            projection=[Hit("projection", "p1", 5, date="2026-03-01"),
                        Hit("projection", "p2", 5, date="2026-02-01"),
                        Hit("projection", "p3", 5, date="2026-01-01")],
            collection=[Hit("collection", "c1", 3), Hit("collection", "c2", 2)],
            genelist=[Hit("genelist", "g1", 1, date="2022-01-01")],
        )
        response = run_search(Criteria(q="hair"), None, None, searchers)
        assert keys(response) == ["p1", "g1", "c1", "p2", "p3", "c2"]

    def test_relevance_without_scores_falls_back_to_newest_first(self):
        searchers = stub_searchers(
            dataset=[Hit("dataset", "old", date="2025-01-01"), Hit("dataset", "undated")],
            study=[Hit("study", "new", date="2026-06-01")],
        )
        response = run_search(Criteria(), None, None, searchers)
        assert keys(response) == ["new", "old", "undated"]

    def test_title_sort_ignores_case(self):
        searchers = stub_searchers(
            dataset=[Hit("dataset", "b", title="beta"), Hit("dataset", "a", title="Alpha")],
            collection=[Hit("collection", "c", title="Gamma")],
        )
        response = run_search(Criteria(sort="title"), None, None, searchers)
        assert keys(response) == ["a", "b", "c"]

    def test_pages_and_hydrates_only_the_current_page(self):
        hits = [Hit("dataset", f"d{i}", date=f"2026-01-0{9 - i}") for i in range(1, 6)]
        searchers = stub_searchers(dataset=hits)
        response = run_search(Criteria(sort="date", page=2, limit=2), None, None, searchers)
        assert keys(response) == ["d3", "d4"]
        assert searchers["dataset"].hydrated == [["d3", "d4"]]
        assert response["pagination"] == {"total_results": 5, "current_page": 2, "limit": 2,
                                          "total_pages": 3, "next_page": 3, "prev_page": 1}

    def test_types_not_on_the_page_are_not_hydrated(self):
        searchers = stub_searchers(dataset=[Hit("dataset", "d1")], genelist=[Hit("genelist", "g1")])
        run_search(Criteria(types=["dataset"]), None, None, searchers)
        assert searchers["genelist"].hydrated == []

    def test_result_that_vanished_before_hydration_is_dropped(self):
        searchers = stub_searchers(dataset=[Hit("dataset", "d1", date="2026-02-01"),
                                            Hit("dataset", "gone", date="2026-01-01")])
        searchers["dataset"].missing = {"gone"}
        response = run_search(Criteria(), None, None, searchers)
        assert keys(response) == ["d1"]

    def test_page_past_the_end_is_empty(self):
        searchers = stub_searchers(dataset=[Hit("dataset", "d1")])
        response = run_search(Criteria(page=4), None, None, searchers)
        assert response["results"] == []
        assert response["pagination"]["total_results"] == 1
        assert response["pagination"]["next_page"] is None

    def test_no_matches(self):
        response = run_search(Criteria(q="nothing"), None, None, stub_searchers())
        assert response["results"] == []
        assert response["counts"]["all"] == 0
        assert response["pagination"]["total_pages"] == 0
