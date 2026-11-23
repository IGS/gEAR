"""
Saved projections in the unified search (#1383).

Each test writes a small projections tree (projections.json plus result CSVs) to a temporary
directory. The dataset and gene list lookups are replaced with stubs that return fixed rows.
"""

import json
import os
import time
import types
from datetime import date, timedelta

import pytest
from gear.search import Criteria
from gear.search import projections
from gear.search.projections import ProjectionSearcher, SavedRun, read_saved_runs

OWNER = types.SimpleNamespace(id=7)
TODAY = date.today().isoformat()


def write_registry(root, dataset_id, registry, csvs=(), raw=None):
    """Write by_dataset/<dataset_id>/projections.json and an empty <uuid>.csv for each uuid in csvs."""
    folder = root / "by_dataset" / dataset_id
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "projections.json").write_text(raw if raw is not None else json.dumps(registry))
    for uuid in csvs:
        (folder / f"{uuid}.csv").write_text("")
    return folder


def set_mtime(path, days_ago):
    stamp = time.time() - days_ago * 86400
    os.utime(path, (stamp, stamp))


def run(uuid, algorithm="pca", zscore=False, day=TODAY):
    return SavedRun(uuid, algorithm, zscore, day)


# ---------------------------------------------------------------- reading the registry


class TestReadSavedRuns:
    def test_runs_are_grouped_by_dataset_and_gene_list(self, tmp_path):
        write_registry(tmp_path, "DS1", {
            "cartA": [{"uuid": "u1", "algorithm": "pca", "zscore": False},
                      {"uuid": "u2", "algorithm": "nmf", "zscore": True}],
            "cartB": [{"uuid": "u3", "algorithm": "fixednmf", "zscore": False}],
        }, csvs=["u1", "u2", "u3"])
        write_registry(tmp_path, "DS2", {"cartA": [{"uuid": "u4", "algorithm": "binary", "zscore": False}]},
                       csvs=["u4"])
        assert read_saved_runs(tmp_path) == {
            ("DS1", "cartA"): [run("u1"), run("u2", "nmf", True)],
            ("DS1", "cartB"): [run("u3", "fixednmf")],
            ("DS2", "cartA"): [run("u4", "binary")],
        }

    def test_run_without_a_result_file_is_left_out(self, tmp_path):
        write_registry(tmp_path, "DS1", {"cartA": [{"uuid": "u1", "algorithm": "pca", "zscore": False},
                                                   {"uuid": "missing", "algorithm": "nmf", "zscore": False}],
                                         "cartB": [{"uuid": "gone", "algorithm": "pca", "zscore": False}]},
                       csvs=["u1"])
        assert read_saved_runs(tmp_path) == {("DS1", "cartA"): [run("u1")]}

    def test_legacy_entries_are_read(self, tmp_path):
        # Older registries key gene lists as "cart.<share_id>", say is_pca instead of algorithm,
        #  and have no zscore
        write_registry(tmp_path, "DS1", {"cart.cartA": [{"uuid": "u1", "is_pca": 1}, {"uuid": "u2", "is_pca": 0}]},
                       csvs=["u1", "u2"])
        assert read_saved_runs(tmp_path) == {("DS1", "cartA"): [run("u1", "pca"), run("u2", "nmf")]}

    def test_repeated_settings_appear_once(self, tmp_path):
        write_registry(tmp_path, "DS1", {
            "cart.cartA": [{"uuid": "u1", "is_pca": 1}],
            "cartA": [{"uuid": "u1", "algorithm": "pca", "zscore": False},
                      {"uuid": "u1", "algorithm": "pca", "zscore": False}],
        }, csvs=["u1"])
        assert read_saved_runs(tmp_path) == {("DS1", "cartA"): [run("u1")]}

    def test_run_date_is_when_its_result_file_was_written(self, tmp_path):
        folder = write_registry(tmp_path, "DS1", {"cartA": [{"uuid": "u1", "algorithm": "pca", "zscore": False}]},
                                csvs=["u1"])
        set_mtime(folder / "u1.csv", days_ago=10)
        expected = (date.today() - timedelta(days=10)).isoformat()
        assert read_saved_runs(tmp_path)[("DS1", "cartA")][0].date == expected

    @pytest.mark.parametrize("raw", ["{}", "{not json", "[]"])
    def test_empty_or_unreadable_registry_is_skipped(self, tmp_path, raw):
        write_registry(tmp_path, "DS1", None, raw=raw)
        write_registry(tmp_path, "DS2", {"cartA": [{"uuid": "u1", "algorithm": "pca", "zscore": False}]},
                       csvs=["u1"])
        assert list(read_saved_runs(tmp_path)) == [("DS2", "cartA")]

    def test_no_projections_directory(self, tmp_path):
        assert read_saved_runs(tmp_path / "nowhere") == {}


# ---------------------------------------------------------------- the searcher

DATASETS = {
    "DS1": {"id": "DS1", "share_id": "ds1share", "title": "Neonatal mouse cochlea", "ldesc": "Hair cell lineage.",
            "organism_id": 1, "organism": "Mouse", "owner_id": 7, "owner_name": "Ann Owner", "is_public": 1,
            "date_added": "2026-01-01", "access": "Yours"},
    "DS2": {"id": "DS2", "share_id": "ds2share", "title": "Adult utricle", "ldesc": "Regeneration time course.",
            "organism_id": 1, "organism": "Mouse", "owner_id": 5, "owner_name": "Bo Other", "is_public": 1,
            "date_added": "2026-01-01", "access": "Public"},
}
CARTS = {
    "cartA": {"share_id": "cartA", "label": "Hair cell markers", "ldesc": "Genes enriched in hair cells."},
    "cartB": {"share_id": "cartB", "label": "Cortex layers", "ldesc": "Layer signatures."},
}


@pytest.fixture
def lookups(monkeypatch):
    """Stub lookups. Set .hidden_datasets / .hidden_carts to hide rows; .calls records each call."""
    stub = types.SimpleNamespace(hidden_datasets=set(), hidden_carts=set(), calls=[])

    def visible_datasets(cursor, user, ids=(), share_ids=(), ownership=(), organism_ids=()):
        stub.calls.append(("datasets", sorted(ids), list(ownership), list(organism_ids)))
        return [DATASETS[i] for i in ids if i in DATASETS and i not in stub.hidden_datasets]

    def visible_gene_carts(cursor, user, share_ids):
        stub.calls.append(("carts", sorted(share_ids)))
        return {s: CARTS[s] for s in share_ids if s in CARTS and s not in stub.hidden_carts}

    monkeypatch.setattr(projections, "visible_datasets", visible_datasets)
    monkeypatch.setattr(projections, "visible_gene_carts", visible_gene_carts)
    return stub


@pytest.fixture
def searcher(tmp_path):
    write_registry(tmp_path, "DS1", {"cartA": [{"uuid": "u1", "algorithm": "pca", "zscore": False},
                                               {"uuid": "u2", "algorithm": "nmf", "zscore": True}],
                                     "cartB": [{"uuid": "u3", "algorithm": "pca", "zscore": False}]},
                   csvs=["u1", "u2", "u3"])
    write_registry(tmp_path, "DS2", {"cartA": [{"uuid": "u4", "algorithm": "pca", "zscore": False}]}, csvs=["u4"])
    return ProjectionSearcher(tmp_path)


def found(searcher, criteria=None, user=None):
    hits = searcher.find(criteria or Criteria(), user, None)
    assert all(hit.type == "projection" for hit in hits)
    return {hit.key: hit for hit in hits}


class TestFind:
    def test_no_keywords_lists_every_visible_pair(self, searcher, lookups):
        hits = found(searcher)
        assert set(hits) == {"DS1/cartA", "DS1/cartB", "DS2/cartA"}
        assert hits["DS1/cartA"].title == "Neonatal mouse cochlea"
        assert hits["DS1/cartA"].date == TODAY

    def test_pair_needs_both_dataset_and_gene_list_visible(self, searcher, lookups):
        lookups.hidden_datasets = {"DS2"}
        lookups.hidden_carts = {"cartB"}
        assert set(found(searcher)) == {"DS1/cartA"}

    def test_dataset_keyword_finds_every_pattern_projected_onto_it(self, searcher, lookups):
        assert set(found(searcher, Criteria(q="cochlea"))) == {"DS1/cartA", "DS1/cartB"}

    def test_gene_list_keyword_finds_every_dataset_it_was_projected_onto(self, searcher, lookups):
        assert set(found(searcher, Criteria(q="markers"))) == {"DS1/cartA", "DS2/cartA"}

    @pytest.mark.parametrize("q", ["projection", "projectR", "patterns"])
    def test_concept_and_tool_names_match_every_pair(self, searcher, lookups, q):
        assert len(found(searcher, Criteria(q=q))) == 3

    def test_title_outranks_description(self, searcher, lookups):
        hits = found(searcher, Criteria(q="utricle regeneration"))
        title_and_desc = hits["DS2/cartA"].score
        hits = found(searcher, Criteria(q="regeneration"))
        assert title_and_desc > hits["DS2/cartA"].score > 0

    def test_filters_are_passed_to_the_dataset_lookup(self, searcher, lookups):
        found(searcher, Criteria(ownership=["yours"], organism_ids=[1]), OWNER)
        assert ("datasets", ["DS1", "DS2"], ["yours"], [1]) in lookups.calls

    def test_date_added_filter_uses_the_newest_run(self, tmp_path, lookups):
        folder = write_registry(tmp_path, "DS1", {"cartA": [{"uuid": "old", "algorithm": "pca", "zscore": False},
                                                            {"uuid": "new", "algorithm": "nmf", "zscore": False}],
                                                  "cartB": [{"uuid": "older", "algorithm": "pca", "zscore": False}]},
                                csvs=["old", "new", "older"])
        set_mtime(folder / "old.csv", days_ago=400)
        set_mtime(folder / "new.csv", days_ago=3)
        set_mtime(folder / "older.csv", days_ago=500)
        hits = found(ProjectionSearcher(tmp_path), Criteria(date_added="month"))
        assert set(hits) == {"DS1/cartA"}
        assert hits["DS1/cartA"].date == (date.today() - timedelta(days=3)).isoformat()

    def test_no_saved_projections_needs_no_database(self, tmp_path, lookups):
        assert found(ProjectionSearcher(tmp_path)) == {}
        assert lookups.calls == []


class TestHydrate:
    def test_result_fields(self, searcher, lookups):
        result = searcher.hydrate(["DS1/cartA"], OWNER, None)["DS1/cartA"]
        # Saved runs open through the projection permalink, as the Projection Tool's own share links do
        first_run = "./p?p=p&s=ds1share&c=cartA&algo=pca&zscore=0"

        assert result["type"] == "projection"
        assert result["id"] == "DS1/cartA"
        assert result["share_id"] == "ds1share"
        assert result["title"] == "Neonatal mouse cochlea"
        assert result["desc"] == "Hair cell lineage."
        assert result["access"] == "Yours"
        assert result["organism"] == "Mouse"
        assert result["date"] == TODAY
        assert result["subline"] == "Projection · Hair cell markers · 2 saved runs"
        assert ["Pattern source", "Hair cell markers"] in result["meta"]
        assert ["Owner", "Ann Owner"] in result["meta"]
        assert result["actions"]["primary"] == {"label": "View in Projection Tool", "url": first_run}
        # The result is the dataset and pattern pair; an action for the pattern alone would be about something else
        assert "secondary" not in result["actions"]
        assert result["actions"]["share_url"] == first_run
        assert result["extra"]["runs"] == [
            {"label": "PCA", "url": first_run},
            {"label": "NMF · z-scored",
             "url": "./p?p=p&s=ds1share&c=cartA&algo=nmf&zscore=1"},
        ]

    def test_one_run_is_singular(self, searcher, lookups):
        assert searcher.hydrate(["DS2/cartA"], None, None)["DS2/cartA"]["subline"] == \
            "Projection · Hair cell markers · 1 saved run"

    def test_tool_name_never_appears(self, searcher, lookups):
        result = searcher.hydrate(["DS1/cartA"], OWNER, None)["DS1/cartA"]
        assert "projectr" not in json.dumps(result).lower()

    def test_pairs_no_longer_visible_or_unknown_are_left_out(self, searcher, lookups):
        lookups.hidden_carts = {"cartB"}
        assert set(searcher.hydrate(["DS1/cartA", "DS1/cartB", "DS9/cartA", "bad-key"], OWNER, None)) == {"DS1/cartA"}
