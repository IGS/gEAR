"""
Collections (layouts) in the unified search (#1383), run against a cursor that records each
query and returns rows chosen by the table the query reads.
"""

import types

import pytest
from gear.search import Criteria
from gear.search.access import DATASET_LIVE_CLAUSE
from gear.search.layouts import DATASET_PREVIEW_COUNT, CollectionSearcher
from gear.utils.fulltext import to_boolean_mode_query, to_like_patterns
from helpers.scripted_cursor import ScriptedCursor

OWNER = types.SimpleNamespace(id=7)

DATASET_MATCH = "MATCH(d.title, d.ldesc, d.geo_id, d.pubmed_id) AGAINST(%s IN BOOLEAN MODE)"


def find(rows=(), criteria=None, user=OWNER):
    cursor = ScriptedCursor({"FROM layout l": list(rows)})
    hits = CollectionSearcher().find(criteria or Criteria(), user, cursor)
    assert all(hit.type == "collection" for hit in hits)
    return {hit.key: hit for hit in hits}, cursor


def only_query(cursor):
    assert len(cursor.queries) == 1
    return cursor.queries[0]


class TestFind:
    def test_no_keywords_lists_every_visible_collection(self):
        # share_id, label, number of the user's visible datasets in it that match the keywords
        hits, cursor = find([("lay1", "Cochlea atlas", 0), ("lay2", "Brain", 0)])
        assert set(hits) == {"lay1", "lay2"}
        assert (hits["lay1"].title, hits["lay1"].score, hits["lay1"].date) == ("Cochlea atlas", 0, None)
        qry, _ = only_query(cursor)
        assert "MATCH(" not in qry and "HAVING" not in qry

    def test_default_visibility(self):
        _, cursor = find()
        qry, params = only_query(cursor)
        assert "l.is_public = 1 OR l.is_domain = 1" in qry and "l.user_id = %s" in qry
        assert "layout_group_membership" in qry
        assert 7 in params

    def test_logged_out_sees_public_and_site_collections_only(self):
        _, cursor = find(criteria=Criteria(ownership=["yours"]), user=None)
        qry, params = only_query(cursor)
        assert "l.is_public = 1 OR l.is_domain = 1" in qry and "l.user_id = %s" not in qry
        assert params == []

    def test_ownership_filter(self):
        _, cursor = find(criteria=Criteria(ownership=["group"]))
        qry, _ = only_query(cursor)
        assert "layout_group_membership" in qry and "l.is_public" not in qry

    def test_keywords_search_the_label_and_the_member_datasets(self):
        _, cursor = find(criteria=Criteria(q="hair potential"))
        qry, params = only_query(cursor)
        assert DATASET_MATCH in qry
        assert qry.count("l.label LIKE %s") == 2 and qry.count("d.title LIKE %s") == 2
        assert "HAVING" in qry
        assert to_boolean_mode_query("hair potential", prefix=True, boost_phrase=True) in params
        assert set(to_like_patterns("hair potential")) <= set(params)

    def test_only_live_datasets_the_user_may_see_are_searched(self):
        # Otherwise a private dataset's title would reveal which collections hold it
        _, cursor = find(criteria=Criteria(q="hair"))
        qry, _ = only_query(cursor)
        assert DATASET_LIVE_CLAUSE in qry
        assert "d.is_public = 1" in qry and "d.owner_id = %s" in qry and "dataset_shares" in qry

    def test_label_match_outranks_a_dataset_match(self):
        hits, _ = find([("label", "Hair cell atlas", 0), ("member", "Inner ear", 1)], Criteria(q="hair"))
        assert hits["label"].score > hits["member"].score > 0

    def test_more_matching_datasets_rank_higher(self):
        hits, _ = find([("one", "Inner ear", 1), ("three", "Organ of Corti", 3)], Criteria(q="hair"))
        assert hits["three"].score > hits["one"].score

    def test_label_match_inside_an_underscore_name_still_scores(self):
        hits, _ = find([("lay1", "Groves_potential_targets", 0)], Criteria(q="potential"))
        assert hits["lay1"].score > 0

    def test_quoted_word_matches_only_whole_words_in_the_label(self):
        # SQL's label check is rough ("%hair%" also finds "Chairs"); the exact rules drop the extra rows
        hits, cursor = find([("whole", "Hair atlas", 0), ("inside", "Chairs", 0)], Criteria(q='"hair"'))
        assert set(hits) == {"whole"}
        assert "l.label LIKE %s" in only_query(cursor)[0]

    def test_organism_filter_needs_a_visible_dataset_of_that_organism(self):
        _, cursor = find(criteria=Criteria(organism_ids=[1, 2]))
        qry, params = only_query(cursor)
        assert "d.organism_id IN (%s, %s)" in qry and "HAVING" in qry
        assert {1, 2} <= set(params)

    def test_date_added_filter_leaves_everything_out(self):
        hits, cursor = find([("lay1", "Cochlea atlas", 0)], Criteria(date_added="week"))
        assert hits == {} and cursor.queries == []


# ---------------------------------------------------------------- hydrate

# id, share_id, label, user_id, user_name, is_public, is_domain
LAYOUT = (21, "lay1", "Cochlea atlas", 7, "Ann Owner", 0, 0)
# collection share_id, dataset share_id, dataset title, organism
MEMBERS = [("lay1", "ds1", "Neonatal cochlea", "Mouse"),
           ("lay1", "ds1", "Neonatal cochlea", "Mouse"),       # a second display of the same dataset
           ("lay1", "ds2", "Human inner ear", "Human"),
           ("lay1", "ds3", "Adult utricle", "Mouse")]


def hydrate(keys, layouts, members=(), user=OWNER):
    cursor = ScriptedCursor({"layout_displays": list(members), "FROM layout l": list(layouts)})
    return CollectionSearcher().hydrate(keys, user, cursor), cursor


class TestHydrate:
    def test_result_fields(self):
        results, _ = hydrate(["lay1"], [LAYOUT], MEMBERS)
        result = results["lay1"]
        url = "./p?l=lay1"
        assert result["type"] == "collection"
        assert result["id"] == result["share_id"] == "lay1"
        assert result["title"] == "Cochlea atlas"
        assert result["desc"] is None
        assert result["access"] == "Yours"
        assert result["organism"] == "Human, Mouse"
        assert result["date"] is None
        assert result["subline"] == "Collection · 3 datasets · Human, Mouse"
        assert ["Owner", "Ann Owner"] in result["meta"]
        assert result["chips"] == ["Human", "Mouse"]
        assert result["extra"] == {
            "datasets": [{"title": "Neonatal cochlea", "url": "./p?s=ds1"},
                         {"title": "Human inner ear", "url": "./p?s=ds2"},
                         {"title": "Adult utricle", "url": "./p?s=ds3"}],
            "dataset_count": 3,
        }
        assert result["actions"] == {"primary": {"label": "Open in Gene Expression", "url": url},
                                     "secondary": None, "share_url": url}

    @pytest.mark.parametrize("layout, user, access", [
        ((21, "lay1", "Mine", 7, "Ann Owner", 0, 0), OWNER, "Yours"),
        ((22, "lay2", "Public", 5, "Bo Other", 1, 0), OWNER, "Public"),
        ((23, "lay3", "Site", 5, "Bo Other", 0, 1), None, "Public"),            # highlighted by the site
        ((24, "lay4", "Group", 5, "Bo Other", 0, 0), OWNER, "Group-affiliated"),
    ])
    def test_access_label(self, layout, user, access):
        results, _ = hydrate([layout[1]], [layout], user=user)
        assert results[layout[1]]["access"] == access

    def test_preview_is_capped_but_count_is_not(self):
        members = [("lay1", f"ds{i}", f"Dataset {i}", "Mouse") for i in range(DATASET_PREVIEW_COUNT + 4)]
        results, _ = hydrate(["lay1"], [LAYOUT], members)
        assert len(results["lay1"]["extra"]["datasets"]) == DATASET_PREVIEW_COUNT
        assert results["lay1"]["extra"]["dataset_count"] == DATASET_PREVIEW_COUNT + 4

    def test_collection_with_no_visible_datasets(self):
        results, _ = hydrate(["lay1"], [LAYOUT])
        result = results["lay1"]
        assert result["subline"] == "Collection · 0 datasets"
        assert result["organism"] is None and result["chips"] == []

    def test_one_dataset_is_singular(self):
        results, _ = hydrate(["lay1"], [LAYOUT], MEMBERS[:1])
        assert results["lay1"]["subline"] == "Collection · 1 dataset · Mouse"

    def test_member_datasets_are_live_and_visible_and_fetched_in_one_query(self):
        _, cursor = hydrate(["lay1", "lay2"], [LAYOUT, (22, "lay2", "Other", 7, "Ann Owner", 0, 0)], MEMBERS,
                            user=None)
        [(qry, params)] = cursor.queries_with("layout_displays")
        assert DATASET_LIVE_CLAUSE in qry and "d.is_public = 1" in qry and "d.owner_id = %s" not in qry
        assert params == ["lay1", "lay2"]

    def test_lookup_uses_default_visibility(self):
        _, cursor = hydrate(["lay1"], [LAYOUT], user=None)
        qry, params = cursor.queries[0]
        assert "l.share_id IN (%s)" in qry and "l.is_public = 1 OR l.is_domain = 1" in qry
        assert params == ["lay1"]

    def test_collections_not_returned_by_the_lookup_are_left_out(self):
        results, _ = hydrate(["lay1", "gone"], [LAYOUT], MEMBERS)
        assert set(results) == {"lay1"}

    def test_nothing_to_hydrate_runs_no_query(self):
        results, cursor = hydrate([], [])
        assert results == {} and cursor.queries == []
