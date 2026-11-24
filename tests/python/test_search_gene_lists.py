"""
Gene lists in the unified search (#1383), run against a cursor that records each query and
returns rows chosen by the table the query reads.
"""

import types
from datetime import datetime

import pytest
from gear.search import Criteria
from gear.search.gene_lists import GENE_PREVIEW_COUNT, GeneListSearcher
from gear.utils.fulltext import to_boolean_mode_query, to_like_patterns
from helpers.scripted_cursor import ScriptedCursor

OWNER = types.SimpleNamespace(id=7)


def find(cursor, criteria=None, user=OWNER, carts_dir="."):
    hits = GeneListSearcher(carts_dir).find(criteria or Criteria(), user, cursor)
    assert all(hit.type == "genelist" for hit in hits)
    return hits


class TestFind:
    def test_rows_become_hits(self):
        cursor = ScriptedCursor({"FROM gene_cart gc": [("cartA", "Hair cell markers", datetime(2026, 3, 4, 5, 6), 2.5)]})
        [hit] = find(cursor, Criteria(q="hair"))
        assert (hit.key, hit.title, hit.date, hit.score) == ("cartA", "Hair cell markers", "2026-03-04", 2.5)

    def test_keywords_use_the_same_full_text_and_label_search_as_the_gene_list_page(self):
        cursor = ScriptedCursor()
        find(cursor, Criteria(q="hair potential"))
        [(qry, params)] = cursor.queries
        assert "MATCH(gc.label, gc.ldesc) AGAINST(%s IN BOOLEAN MODE)" in qry
        assert qry.count("gc.label LIKE %s") == 2
        assert to_boolean_mode_query("hair potential", prefix=True, boost_phrase=True) in params
        assert set(to_like_patterns("hair potential")) <= set(params)

    def test_no_keywords_lists_everything_visible(self):
        cursor = ScriptedCursor({"FROM gene_cart gc": [("cartA", "Hair cell markers", None, 0)]})
        [hit] = find(cursor)
        assert hit.score == 0 and hit.date is None
        assert "MATCH(" not in cursor.queries[0][0]

    def test_label_only_match_still_scores(self):
        # A LIKE match inside an underscore name has no full-text score, but still ranks above nothing
        cursor = ScriptedCursor({"FROM gene_cart gc": [("cartA", "Groves_potential_targets", None, 0)]})
        [hit] = find(cursor, Criteria(q="potential"))
        assert hit.score > 0

    def test_default_visibility(self):
        cursor = ScriptedCursor()
        find(cursor)
        qry, params = cursor.queries[0]
        assert "gc.is_public = 1" in qry and "gc.user_id = %s" in qry and "user_group_membership" in qry
        assert 7 in params

    def test_logged_out_sees_public_only(self):
        cursor = ScriptedCursor()
        find(cursor, Criteria(ownership=["yours"]), user=None)
        qry, params = cursor.queries[0]
        assert "gc.is_public = 1" in qry and "gc.user_id = %s" not in qry
        assert params == []

    def test_ownership_filter(self):
        cursor = ScriptedCursor()
        find(cursor, Criteria(ownership=["yours"]))
        qry, _ = cursor.queries[0]
        assert "gc.user_id = %s" in qry and "gc.is_public = 1" not in qry

    def test_shared_only_matches_nothing(self):
        # Gene lists can't be shared with a user, so "shared" alone has no rule
        cursor = ScriptedCursor()
        find(cursor, Criteria(ownership=["shared"]))
        assert "(0 = 1)" in cursor.queries[0][0]

    def test_organism_filter(self):
        cursor = ScriptedCursor()
        find(cursor, Criteria(organism_ids=[1, 2]))
        qry, params = cursor.queries[0]
        assert "gc.organism_id IN (%s, %s)" in qry
        assert params[-2:] == [1, 2]

    @pytest.mark.parametrize("date_added, interval", [("week", "INTERVAL 1 WEEK"), ("month", "INTERVAL 1 MONTH"),
                                                      ("year", "INTERVAL 1 YEAR")])
    def test_date_added_filter(self, date_added, interval):
        cursor = ScriptedCursor()
        find(cursor, Criteria(date_added=date_added))
        assert f"gc.date_added >= DATE_SUB(NOW(), {interval})" in cursor.queries[0][0]


# ---------------------------------------------------------------- hydrate

# id, share_id, label, ldesc, gctype, organism, user_id, user_name, is_public, date_added
UNWEIGHTED = (11, "cartA", "Hair cell markers", "Genes enriched in hair cells.", "unweighted-list", "Mouse", 7,
              "Ann Owner", 0, datetime(2026, 1, 2, 3, 4))
WEIGHTED = (12, "cartW", "Cochlea patterns", "NMF weights.", "weighted-list", "Human", 5, "Bo Other", 1,
            datetime(2025, 6, 7))
MEMBERS = [(11, "Myo7a"), (11, "Pou4f3"), (11, "Atoh1")]


@pytest.fixture
def carts_dir(tmp_path):
    (tmp_path / "cart.cartW.tab").write_text("identifier\tgene_symbol\tPC1\nENSG1\tGCLC\t0.1\nENSG2\tSOX2\t0.2\n")
    return tmp_path


def hydrate(keys, rows, members=(), user=OWNER, carts_dir="."):
    cursor = ScriptedCursor({"FROM gene_cart_member": list(members), "FROM gene_cart gc": list(rows)})
    return GeneListSearcher(carts_dir).hydrate(keys, user, cursor), cursor


class TestHydrate:
    def test_unweighted_result_fields(self):
        results, _ = hydrate(["cartA"], [UNWEIGHTED], MEMBERS)
        result = results["cartA"]
        assert result["type"] == "genelist"
        assert result["id"] == result["share_id"] == "cartA"
        assert result["title"] == "Hair cell markers"
        assert result["desc"] == "Genes enriched in hair cells."
        assert result["access"] == "Yours"
        assert result["organism"] == "Mouse"
        assert result["date"] == "2026-01-02"
        assert result["subline"] == "Gene list · Mouse · 3 genes"
        assert ["Owner", "Ann Owner"] in result["meta"]
        # Gene symbols keep their stored casing
        assert result["extra"] == {"genes": ["Myo7a", "Pou4f3", "Atoh1"], "gene_count": 3}
        assert result["actions"]["primary"] == {"label": "View in Gene Expression", "url": "./p?c=cartA"}
        assert result["actions"]["secondary"] == {"label": "Open in Gene List Manager",
                                                  "url": "./gene_list_manager.html?share_id=cartA"}
        assert result["actions"]["share_url"] == "./p?c=cartA"

    def test_weighted_genes_come_from_the_cart_file(self, carts_dir):
        results, _ = hydrate(["cartW"], [WEIGHTED], carts_dir=carts_dir)
        result = results["cartW"]
        assert result["subline"] == "Weighted gene list · Human · 2 genes"
        assert result["access"] == "Public"
        assert result["extra"] == {"genes": ["GCLC", "SOX2"], "gene_count": 2}

    def test_missing_cart_file_shows_no_genes(self, tmp_path):
        results, _ = hydrate(["cartW"], [WEIGHTED], carts_dir=tmp_path)
        assert results["cartW"]["extra"] == {"genes": [], "gene_count": None}
        assert results["cartW"]["subline"] == "Weighted gene list · Human"

    def test_preview_is_capped_but_count_is_not(self):
        members = [(11, f"Gene{i}") for i in range(GENE_PREVIEW_COUNT + 5)]
        results, _ = hydrate(["cartA"], [UNWEIGHTED], members)
        assert len(results["cartA"]["extra"]["genes"]) == GENE_PREVIEW_COUNT
        assert results["cartA"]["extra"]["gene_count"] == GENE_PREVIEW_COUNT + 5

    def test_one_gene_is_singular(self):
        results, _ = hydrate(["cartA"], [UNWEIGHTED], MEMBERS[:1])
        assert results["cartA"]["subline"] == "Gene list · Mouse · 1 gene"

    def test_genes_for_the_whole_page_come_from_one_query(self):
        other = (13, "cartB", "Other", "", "unweighted-list", "Mouse", 7, "Ann Owner", 0, None)
        _, cursor = hydrate(["cartA", "cartB"], [UNWEIGHTED, other], MEMBERS)
        assert len(cursor.queries_with("FROM gene_cart_member")) == 1

    def test_lookup_uses_default_visibility(self):
        _, cursor = hydrate(["cartA"], [UNWEIGHTED], MEMBERS, user=None)
        [(qry, params)] = cursor.queries_with("FROM gene_cart gc")
        assert "gc.share_id IN (%s)" in qry and "gc.is_public = 1" in qry
        assert params == ["cartA"]

    def test_logged_out_access_label(self):
        results, _ = hydrate(["cartW"], [WEIGHTED], user=None, carts_dir=".")
        assert results["cartW"]["access"] == "Public"

    def test_group_list_access(self):
        group_row = (14, "cartG", "Group list", "", "unweighted-list", "Mouse", 5, "Bo Other", 0, None)
        results, _ = hydrate(["cartG"], [group_row])
        assert results["cartG"]["access"] == "Group-affiliated"

    def test_lists_not_returned_by_the_lookup_are_left_out(self):
        results, _ = hydrate(["cartA", "gone"], [UNWEIGHTED], MEMBERS)
        assert set(results) == {"cartA"}

    def test_nothing_to_hydrate_runs_no_query(self):
        results, cursor = hydrate([], [])
        assert results == {} and cursor.queries == []
