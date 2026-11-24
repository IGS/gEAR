"""
Datasets in the unified search (#1383).

find() runs against a cursor that records each query. hydrate() gets its rows from a stub of
gear.search.lookups.visible_datasets, and its preview images from a temporary folder.
"""

import types
from datetime import datetime

import pytest
from gear.search import Criteria, datasets
from gear.search.access import DATASET_LIVE_CLAUSE
from gear.search.datasets import DatasetSearcher
from gear.search.fulltext import to_boolean_mode_query, to_like_patterns
from helpers.scripted_cursor import ScriptedCursor

OWNER = types.SimpleNamespace(id=7)


def find(cursor, criteria=None, user=OWNER):
    hits = DatasetSearcher().find(criteria or Criteria(), user, cursor)
    assert all(hit.type == "dataset" for hit in hits)
    return hits


class TestFind:
    def test_rows_become_hits_keyed_by_share_id(self):
        cursor = ScriptedCursor({"FROM dataset d": [("ds1share", "P1, mouse, scRNA-seq, cochlear epithelium (Kelley)",
                                                     datetime(2026, 3, 4, 5, 6), 2.5)]})
        [hit] = find(cursor, Criteria(q="cochlear"))
        assert (hit.key, hit.title, hit.date, hit.score) == (
            "ds1share", "P1, mouse, scRNA-seq, cochlear epithelium (Kelley)", "2026-03-04", 2.5)

    def test_keywords_use_the_same_full_text_and_title_search_as_the_dataset_explorer(self):
        cursor = ScriptedCursor()
        find(cursor, Criteria(q="cochlea potential"))
        [(qry, params)] = cursor.queries
        assert "MATCH(d.title, d.ldesc, d.geo_id, d.pubmed_id) AGAINST(%s IN BOOLEAN MODE)" in qry
        assert qry.count("d.title LIKE %s") == 2
        assert to_boolean_mode_query("cochlea potential", prefix=True, boost_phrase=True) in params
        assert set(to_like_patterns("cochlea potential")) <= set(params)

    def test_no_keywords_lists_everything_visible(self):
        cursor = ScriptedCursor({"FROM dataset d": [("ds1share", "Cochlea atlas", None, 0)]})
        [hit] = find(cursor)
        assert hit.score == 0 and hit.date is None
        assert "MATCH(" not in cursor.queries[0][0]

    def test_title_only_match_still_scores(self):
        # A LIKE match inside an underscore name has no full-text score, but still ranks above nothing
        cursor = ScriptedCursor({"FROM dataset d": [("ds1share", "Groves_potential_targets", None, 0)]})
        [hit] = find(cursor, Criteria(q="potential"))
        assert hit.score > 0

    def test_only_live_datasets(self):
        # Datasets marked for removal or not fully loaded never appear
        cursor = ScriptedCursor()
        find(cursor)
        assert DATASET_LIVE_CLAUSE in cursor.queries[0][0]

    def test_default_visibility_is_public_owned_or_shared(self):
        cursor = ScriptedCursor()
        find(cursor)
        qry, params = cursor.queries[0]
        assert "d.is_public = 1" in qry and "d.owner_id = %s" in qry and "dataset_shares" in qry
        assert "user_group_membership" not in qry
        assert 7 in params

    def test_group_ownership_filter(self):
        cursor = ScriptedCursor()
        find(cursor, Criteria(ownership=["group"]))
        qry, _ = cursor.queries[0]
        assert "user_group_membership" in qry and "d.is_public = 1" not in qry

    def test_logged_out_sees_public_only(self):
        cursor = ScriptedCursor()
        find(cursor, Criteria(ownership=["yours"]), user=None)
        qry, params = cursor.queries[0]
        assert "d.is_public = 1" in qry and "d.owner_id = %s" not in qry
        assert 7 not in params

    def test_organism_filter(self):
        cursor = ScriptedCursor()
        find(cursor, Criteria(organism_ids=[1, 2]))
        qry, params = cursor.queries[0]
        assert "d.organism_id IN (%s, %s)" in qry
        assert {1, 2} <= set(params)

    @pytest.mark.parametrize("date_added, interval", [("week", "INTERVAL 1 WEEK"), ("month", "INTERVAL 1 MONTH"),
                                                      ("year", "INTERVAL 1 YEAR")])
    def test_date_filter(self, date_added, interval):
        cursor = ScriptedCursor()
        find(cursor, Criteria(date_added=date_added))
        assert f"d.date_added >= DATE_SUB(NOW(), {interval})" in cursor.queries[0][0]


# ---------------------------------------------------------------- hydrate


def dataset_row(**overrides):
    """A row as gear.search.lookups.visible_datasets returns it."""
    row = {
        "id": "DS1", "share_id": "ds1share", "title": "P1, mouse, scRNA-seq, cochlear epithelium (Kelley)",
        "ldesc": "Cochlear epithelium at P1.", "organism_id": 1, "organism": "Mouse", "owner_id": 7,
        "owner_name": "Ann Owner", "is_public": 1, "date_added": "2026-01-02", "access": "Yours",
        "dtype": "single-cell-rnaseq", "geo_id": "GSE12345", "pubmed_id": "31234567",
    }
    row.update(overrides)
    return row


@pytest.fixture
def lookup(monkeypatch):
    """Replaces visible_datasets with a stub returning the rows in `lookup.rows` and recording its calls."""
    state = types.SimpleNamespace(rows=[dataset_row()], calls=[])

    def visible_datasets(cursor, user, ids=(), share_ids=(), ownership=(), organism_ids=()):
        state.calls.append({"user": user, "ids": list(ids), "share_ids": list(share_ids),
                            "ownership": list(ownership), "organism_ids": list(organism_ids)})
        return [row for row in state.rows if row["share_id"] in share_ids]

    monkeypatch.setattr(datasets, "visible_datasets", visible_datasets)
    return state


def hydrate(keys, image_root, user=OWNER):
    return DatasetSearcher(image_root=image_root).hydrate(keys, user, cursor=None)


class TestHydrate:
    def test_nothing_to_hydrate_looks_nothing_up(self, lookup, tmp_path):
        assert hydrate([], tmp_path) == {}
        assert lookup.calls == []

    def test_looks_up_the_page_once_with_default_visibility(self, lookup, tmp_path):
        hydrate(["ds1share", "gone"], tmp_path)
        assert lookup.calls == [{"user": OWNER, "ids": [], "share_ids": ["ds1share", "gone"],
                                 "ownership": [], "organism_ids": []}]

    def test_datasets_the_user_cannot_see_are_left_out(self, lookup, tmp_path):
        assert set(hydrate(["ds1share", "gone"], tmp_path)) == {"ds1share"}

    def test_result(self, lookup, tmp_path):
        result = hydrate(["ds1share"], tmp_path)["ds1share"]
        assert result["type"] == "dataset"
        assert result["id"] == result["share_id"] == "ds1share"
        assert result["title"] == "P1, mouse, scRNA-seq, cochlear epithelium (Kelley)"
        assert result["desc"] == "Cochlear epithelium at P1."
        assert result["access"] == "Yours"
        assert result["organism"] == "Mouse"
        assert result["date"] == "2026-01-02"
        assert result["subline"] == "Dataset · Single-cell RNA-seq · Mouse"
        assert result["meta"] == [["Type", "Single-cell RNA-seq"], ["Organism", "Mouse"], ["Owner", "Ann Owner"],
                                  ["Added", "Jan 2, 2026"], ["GEO ID", "GSE12345"], ["PubMed ID", "31234567"]]

    def test_actions_open_the_dataset(self, lookup, tmp_path):
        actions = hydrate(["ds1share"], tmp_path)["ds1share"]["actions"]
        assert actions == {
            "primary": {"label": "View in Gene Expression", "url": "./p?s=ds1share"},
            "secondary": {"label": "Open in Dataset Explorer", "url": "./dataset_explorer.html?share_id=ds1share"},
            "share_url": "./p?s=ds1share",
        }

    def test_geo_and_pubmed_links(self, lookup, tmp_path):
        links = hydrate(["ds1share"], tmp_path)["ds1share"]["links"]
        assert links == [
            {"label": "GEO: GSE12345", "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12345"},
            {"label": "PubMed: 31234567", "url": "https://pubmed.ncbi.nlm.nih.gov/31234567/"},
        ]

    def test_no_geo_or_pubmed_id_means_no_links_or_rows(self, lookup, tmp_path):
        lookup.rows = [dataset_row(geo_id=None, pubmed_id="")]
        result = hydrate(["ds1share"], tmp_path)["ds1share"]
        assert result["links"] == []
        assert [name for name, _ in result["meta"]] == ["Type", "Organism", "Owner", "Added"]

    @pytest.mark.parametrize("dtype, label", [
        ("bulk-rnaseq", "Bulk RNA-seq"),
        ("microarray", "Microarray"),
        ("single-cell-rnaseq", "Single-cell RNA-seq"),
        ("scRNA-seq", "Single-cell RNA-seq"),
        ("single-nucleus-rnaseq", "Single-nucleus RNA-seq"),
        ("gosling", "Epigenetic"),
        ("spatial", "Spatial transcriptomics"),
        ("spatial-h5ad", "Spatial transcriptomics"),
        ("svg-expression", "svg-expression"),   # an unlisted type is shown as stored
    ])
    def test_dataset_type_labels(self, lookup, tmp_path, dtype, label):
        lookup.rows = [dataset_row(dtype=dtype)]
        assert hydrate(["ds1share"], tmp_path)["ds1share"]["meta"][0] == ["Type", label]

    def test_no_dataset_type(self, lookup, tmp_path):
        lookup.rows = [dataset_row(dtype=None)]
        result = hydrate(["ds1share"], tmp_path)["ds1share"]
        assert result["subline"] == "Dataset · Mouse"
        assert "Type" not in [name for name, _ in result["meta"]]


class TestPreviewImage:
    """The same choice of image as the Dataset Explorer (search_datasets.cgi)."""

    def image(self, tmp_path):
        return hydrate(["ds1share"], tmp_path)["ds1share"]["extra"]["preview_image_url"]

    def test_default_display_image_first(self, lookup, tmp_path):
        (tmp_path / "DS1.default.png").touch()
        (tmp_path / "DS1.single.default.png").touch()
        assert self.image(tmp_path) == "./img/dataset_previews/DS1.default.png"

    def test_then_the_single_gene_default(self, lookup, tmp_path):
        (tmp_path / "DS1.single.default.png").touch()
        assert self.image(tmp_path) == "./img/dataset_previews/DS1.single.default.png"

    def test_epigenetic_datasets_without_an_image_use_the_gosling_picture(self, lookup, tmp_path):
        lookup.rows = [dataset_row(dtype="gosling")]
        assert self.image(tmp_path) == "./img/dataset_previews/gosling.png"

    def test_no_image_means_no_preview(self, lookup, tmp_path):
        # The explorer shows a "missing" placeholder; the search preview leaves the section out instead
        assert self.image(tmp_path) is None
