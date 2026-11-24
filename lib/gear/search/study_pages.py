"""
Searcher for study pages (#1383): the pages under www/landing/<slug>/.

The folder keeps the name "landing" because it is part of the pages' URLs; everything else
calls them studies. Each folder describes itself in a small JSON file (added on the
alternative-intro-page branch, also read by bin/generate_intro_stats.py):

    {"portals": ["gear"], "featured": true, "date": "2026-01-31", "title": "...",
     "summary": "...", "pmid": "...", "thumbnail": "thumbnail.jpg", "canonical_url": "https://..."}

Only pages whose "portals" include this site's domain_label (or "*") are searched.

Nothing else needs adding to that file:
- Keywords match the title, the summary and the page's own visible text (authors, abstract).
- Linked collections are read from the page's p?l= and layout_id= links.
- Organisms come from the datasets in those collections.

Study pages have no owner, so they count as public. "date" is when the page was added or
announced, and is used for the "Date added" filter and sort.
"""

import html
import json
import re
import sys
from dataclasses import dataclass, field
from datetime import date, timedelta
from pathlib import Path

from gear.search import Hit
from gear.search.access import DATASET_LIVE_CLAUSE, layout_visibility
from gear.search.text import compile_terms, score

LANDING_DIR = Path(__file__).resolve().parents[3] / "www" / "landing"

# The per-page description file. Kept in one place in case it is renamed (e.g. to study.json).
MANIFEST_NAME = "landing.json"

# How much a search term counts for, by where it was found
TITLE_WEIGHT = 5
SUMMARY_WEIGHT = 3
PAGE_TEXT_WEIGHT = 1

# The "Date added" choices, as days back from today (close to MySQL's INTERVAL 1 <unit>)
DATE_RANGE_DAYS = {"week": 7, "month": 30, "year": 365}

# Page parts that aren't visible text: scripts, styles, the <head>, and comments (including
#  the server-side includes for the header and menu)
_NOT_TEXT = re.compile(r"<(script|style|head|template)\b.*?</\1\s*>|<!--.*?-->", re.IGNORECASE | re.DOTALL)
_TAG = re.compile(r"<[^>]+>")
_WHITESPACE = re.compile(r"\s+")

# Collection links on a page: p?l=<share_id> (also as &l= after other parameters) or layout_id=<share_id>
_COLLECTION_LINK = re.compile(r"(?:[?&]l=|[?&]layout_id=)([A-Za-z0-9_-]+)")

# "van der Valk et al. Single-nucleus atlas ..." -> "van der Valk et al."
_SUMMARY_AUTHORS = re.compile(r"^(\S.{0,80}? et al\.)")


@dataclass
class StudyPage:
    slug: str
    title: str
    summary: str
    pmid: str
    date: str | None
    thumbnail: str
    portals: list[str]
    text: str
    collection_share_ids: list[str] = field(default_factory=list)


def _page_text(page_html):
    """The visible text of a page, as one line."""
    text = _TAG.sub(" ", _NOT_TEXT.sub(" ", page_html))
    return _WHITESPACE.sub(" ", html.unescape(text)).strip()


def _collection_links(page_html):
    """Share IDs of the collections a page links to, once each, in page order."""
    return list(dict.fromkeys(_COLLECTION_LINK.findall(page_html)))


def _valid_date(value):
    """value if it is a YYYY-MM-DD date, otherwise None."""
    try:
        return date.fromisoformat(value).isoformat() if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value or "") else None
    except ValueError:
        return None


def _read_page(manifest_path):
    """Build a StudyPage from a folder's description file and index.html. None if unreadable."""
    folder = manifest_path.parent
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(manifest, dict):
            raise ValueError("not a JSON object")
        index_path = folder / "index.html"
        page_html = index_path.read_text(encoding="utf-8") if index_path.exists() else ""
    except (OSError, ValueError) as err:
        print(f"study_pages: skipping {manifest_path}: {err}", file=sys.stderr)
        return None

    return StudyPage(
        slug=folder.name,
        title=str(manifest.get("title") or folder.name).strip(),
        summary=str(manifest.get("summary") or "").strip(),
        pmid=str(manifest.get("pmid") or "").strip(),
        date=_valid_date(manifest.get("date")),
        # Only a file name in the page's own folder
        thumbnail=Path(str(manifest.get("thumbnail") or "")).name,
        portals=[str(portal) for portal in manifest.get("portals") or []],
        text=_page_text(page_html),
        collection_share_ids=_collection_links(page_html),
    )


def _placeholders(values):
    return ", ".join(["%s"] * len(values))


class StudyPageSearcher:
    """find() and hydrate() for study pages. domain_label comes from www/site_domain_prefs.json."""

    def __init__(self, domain_label, landing_dir=LANDING_DIR):
        self.domain_label = domain_label
        self.landing_dir = Path(landing_dir)

    def _pages(self):
        """This portal's study pages, in folder-name order."""
        pages = []
        for manifest_path in sorted(self.landing_dir.glob(f"*/{MANIFEST_NAME}")):
            page = _read_page(manifest_path)
            if page and ("*" in page.portals or self.domain_label in page.portals):
                pages.append(page)
        return pages

    def _visible_collections(self, share_ids, user, cursor):
        """{share_id: label} for the collections among share_ids that the user may see."""
        if not share_ids:
            return {}
        clause, params = layout_visibility(user, [])
        cursor.execute(
            f"SELECT l.share_id, l.label FROM layout l WHERE l.share_id IN ({_placeholders(share_ids)}) AND {clause}",
            [*share_ids, *params],
        )
        return {share_id: label for share_id, label in cursor.fetchall()}

    def _collection_organisms(self, share_ids, user, cursor):
        """{share_id: set of organism ids} of the datasets in the collections the user may see."""
        if not share_ids:
            return {}
        clause, params = layout_visibility(user, [])
        cursor.execute(
            f"""
            SELECT DISTINCT l.share_id, d.organism_id
              FROM layout l
                   JOIN layout_displays ld ON ld.layout_id = l.id
                   JOIN dataset_display dd ON dd.id = ld.display_id
                   JOIN dataset d ON d.id = dd.dataset_id
             WHERE l.share_id IN ({_placeholders(share_ids)}) AND {clause} AND {DATASET_LIVE_CLAUSE}
            """,
            [*share_ids, *params],
        )
        organisms = {}
        for share_id, organism_id in cursor.fetchall():
            organisms.setdefault(share_id, set()).add(organism_id)
        return organisms

    def find(self, criteria, user, cursor):
        # Study pages count as public. Logged-out searches ignore the ownership filter, as they
        #  do for datasets.
        if user and criteria.ownership and "public" not in criteria.ownership:
            return []

        pages = self._pages()

        if criteria.date_added:
            since = (date.today() - timedelta(days=DATE_RANGE_DAYS[criteria.date_added])).isoformat()
            pages = [page for page in pages if page.date and page.date >= since]

        if criteria.organism_ids:
            share_ids = list(dict.fromkeys(sid for page in pages for sid in page.collection_share_ids))
            organisms = self._collection_organisms(share_ids, user, cursor)
            wanted = set(criteria.organism_ids)
            pages = [page for page in pages
                     if any(organisms.get(sid, set()) & wanted for sid in page.collection_share_ids)]

        patterns = compile_terms(criteria.q)
        hits = []
        for page in pages:
            page_score = 0
            if patterns:
                page_score = score(patterns, [(page.title, TITLE_WEIGHT), (page.summary, SUMMARY_WEIGHT),
                                              (page.text, PAGE_TEXT_WEIGHT)])
                if not page_score:
                    continue
            hits.append(Hit("study", page.slug, page_score, page.title, page.date))
        return hits

    def hydrate(self, keys, user, cursor):
        pages = {page.slug: page for page in self._pages()}
        wanted = [pages[key] for key in keys if key in pages]
        if not wanted:
            return {}

        share_ids = list(dict.fromkeys(sid for page in wanted for sid in page.collection_share_ids))
        collection_labels = self._visible_collections(share_ids, user, cursor)
        collection_organisms = self._collection_organisms(list(collection_labels), user, cursor)

        organism_ids = sorted({oid for ids in collection_organisms.values() for oid in ids})
        organism_labels = {}
        if organism_ids:
            cursor.execute(f"SELECT id, label FROM organism WHERE id IN ({_placeholders(organism_ids)})",
                           organism_ids)
            organism_labels = dict(cursor.fetchall())

        return {page.slug: self._result(page, collection_labels, collection_organisms, organism_labels)
                for page in wanted}

    def _result(self, page, collection_labels, collection_organisms, organism_labels):
        page_url = f"./landing/{page.slug}/"
        visible_ids = [sid for sid in page.collection_share_ids if sid in collection_labels]
        collections = [{"label": collection_labels[sid], "url": f"./p?l={sid}"} for sid in visible_ids]

        organism_ids = {oid for sid in visible_ids for oid in collection_organisms.get(sid, ())}
        organisms = sorted(organism_labels[oid] for oid in organism_ids if oid in organism_labels)

        match = _SUMMARY_AUTHORS.match(page.summary)
        authors = match.group(1) if match else None

        collection_count = None
        if collections:
            collection_count = f"{len(collections)} linked collection{'s' if len(collections) != 1 else ''}"

        added = None
        if page.date:
            day = date.fromisoformat(page.date)
            added = f"{day:%b} {day.day}, {day.year}"

        meta = [["Authors", authors], ["Organisms", ", ".join(organisms)],
                ["Linked collections", str(len(collections)) if collections else None], ["Added", added]]

        return {
            "type": "study",
            "id": page.slug,
            "share_id": None,
            "title": page.title,
            "desc": page.summary,
            "access": "Public",
            "organism": ", ".join(organisms) or None,
            "date": page.date,
            "subline": " · ".join(part for part in ("Study", authors, collection_count) if part),
            "meta": [row for row in meta if row[1]],
            "chips": organisms,
            "links": [{"label": "PubMed", "url": f"https://pubmed.ncbi.nlm.nih.gov/{page.pmid}/"}]
                     if page.pmid.isdigit() else [],
            "actions": {
                "primary": {"label": "Open study page", "url": page_url},
                "secondary": {"label": "View linked collection", "url": collections[0]["url"]} if collections else None,
                "share_url": page_url,
            },
            "extra": {
                "thumbnail_url": f"{page_url}{page.thumbnail}" if page.thumbnail else None,
                "collections": collections,
            },
        }
