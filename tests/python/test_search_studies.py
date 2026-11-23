"""
Study pages (the pages under /landing/) in the unified search (#1383).

Each test writes a small www/landing/ tree (landing.json + index.html per page) to a temporary
directory. A fake database cursor answers the collection and organism lookups.
"""

import json
import types
from datetime import date, timedelta

import pytest
from gear.search import Criteria
from gear.search.study_pages import StudyPageSearcher

OWNER = types.SimpleNamespace(id=7)


class FakeCursor:
    """
    Answers the three lookups the searcher makes:
    - organisms of the datasets in some collections (layout_organisms: share_id -> organism ids),
    - labels of the collections the user may see (layouts: share_id -> label),
    - organism labels.
    """

    ORGANISMS = {1: "Mouse", 2: "Human", 3: "Zebrafish", 5: "Chicken"}

    def __init__(self, layouts=None, layout_organisms=None):
        self.layouts = layouts or {}
        self.layout_organisms = layout_organisms or {}
        self.queries = []
        self._rows = []

    def execute(self, qry, params=()):
        params = list(params)
        self.queries.append((qry, params))
        if "dataset_display" in qry:
            self._rows = [(share_id, organism_id) for share_id, ids in self.layout_organisms.items()
                          if share_id in params for organism_id in ids]
        elif "FROM layout" in qry:
            self._rows = [(share_id, label) for share_id, label in self.layouts.items() if share_id in params]
        elif "FROM organism" in qry:
            self._rows = [(id, label) for id, label in self.ORGANISMS.items() if id in params]
        else:
            raise AssertionError(f"unexpected query: {qry}")

    def fetchall(self):
        return self._rows

    def __iter__(self):
        return iter(self._rows)


@pytest.fixture
def landing(tmp_path):
    """add(slug, body=..., **landing_json) writes one study page; .searcher() searches them."""
    class Landing:
        def add(self, slug, body="", **fields):
            folder = tmp_path / slug
            folder.mkdir()
            manifest = {"portals": ["gear"], "featured": True, "date": "2026-01-01", "title": slug,
                        "summary": "", "pmid": "", "thumbnail": "", "canonical_url": ""}
            manifest.update(fields)
            (folder / "landing.json").write_text(json.dumps(manifest))
            (folder / "index.html").write_text(f"<html><body><section id='content-c'>{body}</section></body></html>")

        def searcher(self, domain="gear"):
            return StudyPageSearcher(domain, landing_dir=tmp_path)

    return Landing()


def found(searcher, criteria=None, user=None, cursor=None):
    hits = searcher.find(criteria or Criteria(), user, cursor or FakeCursor())
    assert all(hit.type == "study" for hit in hits)
    return {hit.key: hit for hit in hits}


class TestFind:
    def test_no_keywords_lists_this_portals_pages(self, landing):
        landing.add("ours")
        landing.add("theirs", portals=["nemo"])
        landing.add("everywhere", portals=["*"])
        assert set(found(landing.searcher())) == {"ours", "everywhere"}

    def test_folder_without_landing_json_is_ignored(self, landing, tmp_path):
        landing.add("listed")
        (tmp_path / "draft").mkdir()
        (tmp_path / "draft" / "index.html").write_text("<p>draft</p>")
        assert set(found(landing.searcher())) == {"listed"}

    def test_malformed_landing_json_is_skipped(self, landing, tmp_path):
        landing.add("good")
        landing.add("broken")
        (tmp_path / "broken" / "landing.json").write_text("{not json")
        assert set(found(landing.searcher())) == {"good"}

    def test_words_match_title_summary_and_page_text(self, landing):
        landing.add("t", title="Atlas of the cochlea")
        landing.add("s", summary="Cells of the cochlea.")
        landing.add("b", body="<p id='author-list'>Jane Cochlea, John Smith</p>")
        landing.add("none", title="Neocortex", body="<p>cortex</p>")
        assert set(found(landing.searcher(), Criteria(q="cochlea"))) == {"t", "s", "b"}

    def test_markup_and_scripts_are_not_searched(self, landing):
        landing.add("p", body="<script>const cochlea = 1;</script><a class='cochlea' href='#'>link</a>")
        assert found(landing.searcher(), Criteria(q="cochlea")) == {}

    def test_words_match_as_prefixes(self, landing):
        landing.add("p", title="Cochlear development")
        assert set(found(landing.searcher(), Criteria(q="cochle"))) == {"p"}

    def test_any_word_matches(self, landing):
        landing.add("hair", title="Hair cell atlas")
        landing.add("brain", title="Brain atlas")
        landing.add("other", title="Pancreas")
        assert set(found(landing.searcher(), Criteria(q="hair brain"))) == {"hair", "brain"}

    def test_quoted_text_matches_only_as_a_phrase(self, landing):
        landing.add("phrase", title="Hair cell atlas")
        landing.add("apart", title="Cell types near hair follicles")
        assert set(found(landing.searcher(), Criteria(q='"hair cell"'))) == {"phrase"}

    def test_title_outranks_summary_outranks_page_text(self, landing):
        landing.add("title", title="HIEDRA atlas")
        landing.add("summary", title="Other", summary="Built on HIEDRA.")
        landing.add("body", title="Another", body="<p>HIEDRA</p>")
        hits = found(landing.searcher(), Criteria(q="hiedra"))
        assert hits["title"].score > hits["summary"].score > hits["body"].score > 0

    def test_hit_has_title_and_date(self, landing):
        landing.add("p", title="An atlas", date="2025-08-19")
        hit = found(landing.searcher())["p"]
        assert (hit.title, hit.date) == ("An atlas", "2025-08-19")

    def test_malformed_date_is_treated_as_undated(self, landing):
        landing.add("p", date="August 2025")
        assert found(landing.searcher())["p"].date is None

    def test_date_added_filter_uses_the_page_date(self, landing):
        today = date.today()
        landing.add("recent", date=(today - timedelta(days=3)).isoformat())
        landing.add("old", date=(today - timedelta(days=400)).isoformat())
        landing.add("undated", date="")
        assert set(found(landing.searcher(), Criteria(date_added="month"))) == {"recent"}

    def test_organism_filter_uses_the_linked_collections_datasets(self, landing):
        landing.add("mouse", body="<a href='https://umgear.org/p?l=MouseLayout&g=Sox2'>LINK</a>")
        landing.add("human", body="<a href='./expression.html?layout_id=HumanLayout'>LINK</a>")
        landing.add("no_links")
        cursor = FakeCursor(layout_organisms={"MouseLayout": [1], "HumanLayout": [2]})
        assert set(found(landing.searcher(), Criteria(organism_ids=[1, 3]), cursor=cursor)) == {"mouse"}

    def test_no_organism_filter_needs_no_database(self, landing):
        landing.add("p")
        cursor = FakeCursor()
        found(landing.searcher(), Criteria(q="p"), cursor=cursor)
        assert cursor.queries == []

    @pytest.mark.parametrize("user, ownership, listed", [
        (OWNER, ["public"], True),
        (OWNER, ["yours", "group"], False),
        (OWNER, [], True),
        (None, ["yours"], True),     # logged out, the ownership filter is ignored
    ])
    def test_pages_count_as_public(self, landing, user, ownership, listed):
        landing.add("p")
        hits = found(landing.searcher(), Criteria(ownership=ownership), user)
        assert ("p" in hits) == listed


class TestHydrate:
    BODY = """
        <p id="author-list">WH van der Valk, H Locher</p>
        <a href="https://umgear.org/p?l=L1&g=Sox2">1</a>
        <a href="https://umgear.org/p?l=L2&g=Myo7a">2</a>
        <a href="https://umgear.org/p?l=L1&g=Pou4f3">same collection again</a>
        <a href="https://umgear.org/p?l=Hidden">3</a>
    """

    @pytest.fixture
    def hiedra(self, landing):
        landing.add("hiedra", body=self.BODY, title="HIEDRA atlas", date="2026-02-25",
                    summary="van der Valk et al. Single-nucleus atlas of the developing human inner ear.",
                    pmid="42717251", thumbnail="thumbnail.jpg")
        return landing.searcher()

    @pytest.fixture
    def cursor(self):
        return FakeCursor(layouts={"L1": "HIEDRA overview", "L2": "HIEDRA pseudotime"},
                          layout_organisms={"L1": [2], "L2": [2], "Hidden": [1]})

    def test_result_fields(self, hiedra, cursor):
        result = hiedra.hydrate(["hiedra"], None, cursor)["hiedra"]

        assert result["type"] == "study"
        assert result["id"] == "hiedra"
        assert result["title"] == "HIEDRA atlas"
        assert result["desc"] == "van der Valk et al. Single-nucleus atlas of the developing human inner ear."
        assert result["access"] == "Public"
        assert result["date"] == "2026-02-25"
        assert result["organism"] == "Human"
        assert result["subline"] == "Study · van der Valk et al. · 2 linked collections"
        assert result["actions"]["primary"] == {"label": "Open study page", "url": "./landing/hiedra/"}
        assert result["actions"]["share_url"] == "./landing/hiedra/"
        assert result["links"] == [{"label": "PubMed", "url": "https://pubmed.ncbi.nlm.nih.gov/42717251/"}]
        assert result["extra"]["thumbnail_url"] == "./landing/hiedra/thumbnail.jpg"

    def test_lists_each_visible_linked_collection_once(self, hiedra, cursor):
        result = hiedra.hydrate(["hiedra"], OWNER, cursor)["hiedra"]

        assert result["extra"]["collections"] == [
            {"label": "HIEDRA overview", "url": "./p?l=L1"},
            {"label": "HIEDRA pseudotime", "url": "./p?l=L2"},
        ]
        layout_qry, params = next(q for q in cursor.queries if "FROM layout" in q[0] and "dataset_display" not in q[0])
        assert "l.is_public = 1" in layout_qry      # the collection visibility clause is applied
        assert 7 in params

    def test_organisms_come_only_from_visible_collections(self, hiedra, cursor):
        # "Hidden" holds mouse data, but the user can't see that collection
        assert hiedra.hydrate(["hiedra"], None, cursor)["hiedra"]["organism"] == "Human"

    def test_page_without_pmid_thumbnail_or_collections(self, landing):
        landing.add("hrp", title="Hearing Restoration Project", summary="Datasets and tools from the HRP.")
        result = landing.searcher().hydrate(["hrp"], None, FakeCursor())["hrp"]
        assert result["subline"] == "Study"
        assert result["organism"] is None
        assert result["links"] == []
        assert result["extra"] == {"thumbnail_url": None, "collections": []}

    def test_unknown_key_is_left_out(self, landing):
        landing.add("p")
        assert landing.searcher().hydrate(["gone"], None, FakeCursor()) == {}
