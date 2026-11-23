"""
Cell-cell communication analyses in the unified search (#1383).

Each test writes analysis files (an empty .parquet, optionally with a .json description) to a
temporary directory. The dataset lookup is replaced with a stub that returns fixed rows.
"""

import json
import types

import pytest
from gear.search import Criteria
from gear.search import cell_cell
from gear.search.cell_cell import Analysis, CellCellSearcher, dataset_ref_from_url, read_analyses

OWNER = types.SimpleNamespace(id=7)


def add_analysis(folder, name, info=None):
    (folder / f"{name}.parquet").touch()
    if info is not None:
        (folder / f"{name}.json").write_text(json.dumps(info))
    return f"{name}.parquet"


class TestReadAnalyses:
    def test_lists_parquet_files_with_their_descriptions(self, tmp_path):
        add_analysis(tmp_path, "b_net", {"Title": " Signaling in the cochlea ", "PMID": 123})
        add_analysis(tmp_path, "a_net")
        (tmp_path / "notes.txt").touch()
        analyses = read_analyses(tmp_path)
        assert [a.filename for a in analyses] == ["a_net.parquet", "b_net.parquet"]
        assert analyses[0].info is None
        # Read by gear.dominosignal: trimmed, numbers as strings, blanks as None
        assert analyses[1].info["title"] == "Signaling in the cochlea"
        assert analyses[1].info["pmid"] == "123"
        assert analyses[1].info["author"] is None

    def test_missing_directory(self, tmp_path):
        assert read_analyses(tmp_path / "nowhere") == []


@pytest.mark.parametrize("url, ref", [
    ("https://umgear.org/p?s=abc123", ("share_id", "abc123")),
    ("https://umgear.org/p?share_id=abc123&g=Sox2", ("share_id", "abc123")),
    ("umgear.org/p?s=abc123", ("share_id", "abc123")),
    ("https://umgear.org/expression.html?dataset_id=1f2e3d4c-0000-1111-2222-333344445555",
     ("id", "1f2e3d4c-0000-1111-2222-333344445555")),
    ("https://umgear.org/p?l=layout1", None),     # a collection, not a dataset
    ("umgear.org", None),
    ("", None),
    (None, None),
])
def test_dataset_ref_from_url(url, ref):
    assert dataset_ref_from_url(url) == ref


# ---------------------------------------------------------------- the searcher

DATASET = {"id": "DS1", "share_id": "ds1share", "title": "Adult mouse cochlear epithelium", "ldesc": "",
           "organism_id": 1, "organism": "Mouse", "owner_id": 5, "owner_name": "Bo Other", "is_public": 0,
           "date_added": "2026-01-01", "access": "Group-affiliated"}

LINKED = {"Title": "Signaling between hair cells and supporting cells", "Author": "A Smith", "Year": "2025",
          "DOI": "10.1000/xyz", "PMID": "424242", "gEAR Dataset URL": "https://umgear.org/p?s=ds1share",
          "Description": "Ligand pairs in the organ of Corti."}
UNLINKED = {"Title": "Pancreas niche", "Author": "B Jones", "Year": "2024", "gEAR Dataset URL": "umgear.org",
            "Description": "Stromal signals."}


@pytest.fixture
def lookups(monkeypatch):
    """Stub dataset lookup. Set .visible = False to hide DS1; .calls records each call."""
    stub = types.SimpleNamespace(visible=True, calls=[])

    def visible_datasets(cursor, user, ids=(), share_ids=(), ownership=(), organism_ids=()):
        stub.calls.append({"ids": sorted(ids), "share_ids": sorted(share_ids),
                           "ownership": list(ownership), "organism_ids": list(organism_ids)})
        matches = DATASET["id"] in ids or DATASET["share_id"] in share_ids
        return [DATASET] if stub.visible and matches else []

    monkeypatch.setattr(cell_cell, "visible_datasets", visible_datasets)
    return stub


@pytest.fixture
def searcher(tmp_path):
    add_analysis(tmp_path, "linked_net", LINKED)
    add_analysis(tmp_path, "unlinked_net", UNLINKED)
    add_analysis(tmp_path, "bare_net")
    return CellCellSearcher(tmp_path)


def found(searcher, criteria=None, user=None):
    hits = searcher.find(criteria or Criteria(), user, None)
    assert all(hit.type == "ccc" for hit in hits)
    return {hit.key: hit for hit in hits}


class TestFind:
    def test_no_keywords_lists_every_visible_analysis(self, searcher, lookups):
        hits = found(searcher)
        assert set(hits) == {"linked_net.parquet", "unlinked_net.parquet", "bare_net.parquet"}
        assert hits["linked_net.parquet"].title == "Signaling between hair cells and supporting cells"
        assert all(hit.date is None for hit in hits.values())

    def test_linked_analysis_is_hidden_when_its_dataset_is_not_visible(self, searcher, lookups):
        lookups.visible = False
        assert set(found(searcher)) == {"unlinked_net.parquet", "bare_net.parquet"}

    @pytest.mark.parametrize("q, expected", [
        ("corti", {"linked_net.parquet"}),                   # description
        ("smith", {"linked_net.parquet"}),                   # author
        ("424242", {"linked_net.parquet"}),                  # PMID
        ("epithelium", {"linked_net.parquet"}),              # the linked dataset's title
        ("pancreas", {"unlinked_net.parquet"}),              # title
    ])
    def test_keywords_match_the_description_fields(self, searcher, lookups, q, expected):
        assert set(found(searcher, Criteria(q=q))) == expected

    @pytest.mark.parametrize("q", ["dominoSignal", "ligand", '"cell-cell communication"'])
    def test_concept_and_tool_names_match_every_analysis(self, searcher, lookups, q):
        assert len(found(searcher, Criteria(q=q))) == 3

    def test_title_outranks_description(self, searcher, lookups):
        hits = found(searcher, Criteria(q="stromal pancreas"))
        both = hits["unlinked_net.parquet"].score
        assert both > found(searcher, Criteria(q="stromal"))["unlinked_net.parquet"].score > 0

    def test_filters_are_passed_to_the_dataset_lookup(self, searcher, lookups):
        found(searcher, Criteria(ownership=["group"], organism_ids=[1]), OWNER)
        assert lookups.calls[-1] == {"ids": [], "share_ids": ["ds1share"], "ownership": ["group"],
                                     "organism_ids": [1]}

    def test_unlinked_analyses_count_as_public(self, searcher, lookups):
        assert "unlinked_net.parquet" in found(searcher, Criteria(ownership=["public"]), OWNER)
        assert "unlinked_net.parquet" not in found(searcher, Criteria(ownership=["yours"]), OWNER)

    def test_organism_filter_leaves_out_unlinked_analyses(self, searcher, lookups):
        assert set(found(searcher, Criteria(organism_ids=[1]))) == {"linked_net.parquet"}

    def test_date_added_filter_leaves_everything_out(self, searcher, lookups):
        assert found(searcher, Criteria(date_added="year")) == {}


class TestHydrate:
    def test_linked_result_fields(self, searcher, lookups):
        result = searcher.hydrate(["linked_net.parquet"], OWNER, None)["linked_net.parquet"]
        assert result["type"] == "ccc"
        assert result["id"] == "linked_net.parquet"
        assert result["title"] == "Signaling between hair cells and supporting cells"
        assert result["desc"] == "Ligand pairs in the organ of Corti."
        assert result["access"] == "Group-affiliated"
        assert result["organism"] == "Mouse"
        assert result["subline"] == "Cell-cell communication · Mouse · A Smith, 2025"
        assert ["Dataset", "Adult mouse cochlear epithelium"] in result["meta"]
        assert result["links"] == [{"label": "Read paper", "url": "https://doi.org/10.1000/xyz"},
                                   {"label": "PubMed", "url": "https://pubmed.ncbi.nlm.nih.gov/424242/"}]
        assert result["actions"]["primary"] == {"label": "View signaling results",
                                                "url": "./cell_cell.html?domino_file=linked_net.parquet"}
        assert result["actions"]["secondary"] == {"label": "View in Gene Expression", "url": "./p?s=ds1share"}
        assert result["actions"]["share_url"] == "./cell_cell.html?domino_file=linked_net.parquet"

    def test_unlinked_result_is_public_with_no_dataset_action(self, searcher, lookups):
        result = searcher.hydrate(["unlinked_net.parquet"], None, None)["unlinked_net.parquet"]
        assert result["access"] == "Public"
        assert result["organism"] is None
        assert result["actions"]["secondary"] is None

    def test_file_name_is_url_encoded(self, tmp_path, lookups):
        add_analysis(tmp_path, "net p21&x")
        result = CellCellSearcher(tmp_path).hydrate(["net p21&x.parquet"], None, None)["net p21&x.parquet"]
        assert result["actions"]["primary"]["url"] == "./cell_cell.html?domino_file=net%20p21%26x.parquet"

    def test_tool_name_never_appears_in_visible_text(self, searcher, lookups):
        # The file names contain the tool's name, so check everything except the link URLs
        for key in ("linked_net.parquet", "unlinked_net.parquet", "bare_net.parquet"):
            result = searcher.hydrate([key], OWNER, None)[key]
            visible = [result["title"], result["desc"], result["subline"], result["meta"], result["chips"],
                       [link["label"] for link in result["links"]],
                       [action["label"] for action in (result["actions"]["primary"], result["actions"]["secondary"])
                        if action]]
            assert "domino" not in json.dumps(visible).lower(), key

    def test_untitled_analysis_gets_a_plain_title(self, searcher, lookups):
        result = searcher.hydrate(["bare_net.parquet"], None, None)["bare_net.parquet"]
        assert result["title"] == "Untitled cell-cell communication analysis"

    def test_unknown_or_hidden_analyses_are_left_out(self, searcher, lookups):
        lookups.visible = False
        assert set(searcher.hydrate(["linked_net.parquet", "nope.parquet", "bare_net.parquet"], OWNER, None)) == \
            {"bare_net.parquet"}
