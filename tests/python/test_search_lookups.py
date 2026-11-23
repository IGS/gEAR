"""
Unified search (#1383): the shared dataset and gene list lookups, run against a cursor that
records each query and returns fixed rows.
"""

import types
from datetime import datetime

import pytest
from gear.search import lookups
from gear.search.access import DATASET_LIVE_CLAUSE

OWNER = types.SimpleNamespace(id=7)


class RecordingCursor:
    def __init__(self, rows=()):
        self.rows = list(rows)
        self.queries = []

    def execute(self, qry, params=()):
        self.queries.append((qry, list(params)))

    def fetchall(self):
        return self.rows


def only_query(cursor):
    assert len(cursor.queries) == 1
    qry, params = cursor.queries[0]
    assert qry.count("%s") == len(params)
    return qry, params


@pytest.mark.parametrize("user, owner_id, is_public, is_shared, label", [
    (None, 5, 1, 0, "Public"),
    (OWNER, 7, 1, 0, "Yours"),           # owning it wins over everything else
    (OWNER, 5, 1, 1, "Public"),
    (OWNER, 5, 0, 1, "Shared with you"),
    (OWNER, 5, 0, 0, "Group-affiliated"),  # visible, private, not owned or shared: through a group
])
def test_access_label(user, owner_id, is_public, is_shared, label):
    assert lookups.access_label(user, owner_id, is_public, is_shared) == label


class TestVisibleDatasets:
    ROW = ("DS1", "ds1share", "Cochlea atlas", "Hair cells.", 1, "Mouse", 7, "Ann Owner", 0,
           datetime(2026, 1, 2, 3, 4), "single-cell-rnaseq", "GSE12345", "31234567", 0)

    def test_nothing_to_look_up_runs_no_query(self):
        cursor = RecordingCursor()
        assert lookups.visible_datasets(cursor, OWNER) == []
        assert cursor.queries == []

    def test_query_applies_ids_live_visibility_and_filters(self):
        cursor = RecordingCursor()
        lookups.visible_datasets(cursor, OWNER, ids=["DS1", "DS2"], share_ids=["s3"],
                                 ownership=["yours"], organism_ids=[1, 2])
        qry, params = only_query(cursor)
        assert "d.id IN (%s, %s)" in qry and "d.share_id IN (%s)" in qry
        assert DATASET_LIVE_CLAUSE in qry
        assert "d.owner_id = %s" in qry and "d.is_public = 1" not in qry
        assert "d.organism_id IN (%s, %s)" in qry
        assert {"DS1", "DS2", "s3", 7, 1, 2} <= set(params)

    def test_logged_out_sees_public_only(self):
        cursor = RecordingCursor()
        lookups.visible_datasets(cursor, None, ids=["DS1"], ownership=["yours"])
        qry, params = only_query(cursor)
        assert "d.is_public = 1" in qry and "d.owner_id = %s" not in qry
        assert params == ["DS1"]

    def test_rows_become_dicts_with_an_access_label(self):
        cursor = RecordingCursor(rows=[self.ROW])
        assert lookups.visible_datasets(cursor, OWNER, ids=["DS1"]) == [{
            "id": "DS1", "share_id": "ds1share", "title": "Cochlea atlas", "ldesc": "Hair cells.",
            "organism_id": 1, "organism": "Mouse", "owner_id": 7, "owner_name": "Ann Owner",
            "is_public": 0, "date_added": "2026-01-02", "access": "Yours",
            "dtype": "single-cell-rnaseq", "geo_id": "GSE12345", "pubmed_id": "31234567",
        }]


class TestVisibleGeneCarts:
    def test_nothing_to_look_up_runs_no_query(self):
        cursor = RecordingCursor()
        assert lookups.visible_gene_carts(cursor, OWNER, []) == {}
        assert cursor.queries == []

    def test_query_uses_the_default_gene_list_visibility(self):
        cursor = RecordingCursor(rows=[("cartA", "Hair cell markers", "Markers.")])
        result = lookups.visible_gene_carts(cursor, OWNER, ["cartA", "cartB"])
        qry, params = only_query(cursor)
        assert "gc.share_id IN (%s, %s)" in qry
        assert "gc.is_public = 1" in qry and "gc.user_id = %s" in qry
        assert {"cartA", "cartB", 7} <= set(params)
        assert result == {"cartA": {"share_id": "cartA", "label": "Hair cell markers", "ldesc": "Markers."}}
