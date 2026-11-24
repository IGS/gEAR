"""
Unified search across gEAR entities (#1383, Stage 1).

Each entity type has a searcher module with two functions:

- find(criteria, user, cursor) -> list[Hit]
  Returns every match the user can see: just enough to count, sort and page the results.
- hydrate(keys, user, cursor) -> dict[key, dict]
  Builds the full result for each key, for the current page only.

run_search() combines the searchers: it counts matches per type, applies the type filter,
sorts, pages, and hydrates only the page being returned.
"""

from dataclasses import dataclass, field
from math import ceil

# Order of the type filter buttons in the UI
TYPES = ("dataset", "collection", "genelist", "study", "projection", "ccc")

SORTS = ("relevance", "title", "date")
OWNERSHIPS = ("yours", "shared", "group", "public")
DATE_RANGES = ("week", "month", "year")

DEFAULT_LIMIT = 20
MAX_LIMIT = 100


@dataclass
class Hit:
    """One match, before hydration. date is an ISO "YYYY-MM-DD" string or None."""
    type: str
    key: str
    score: float = 0.0
    title: str = ""
    date: str | None = None


@dataclass
class Criteria:
    q: str = ""
    types: list[str] = field(default_factory=list)       # empty means every type
    organism_ids: list[int] = field(default_factory=list)
    ownership: list[str] = field(default_factory=list)
    date_added: str | None = None
    sort: str = "relevance"
    page: int = 1
    limit: int = DEFAULT_LIMIT


def _csv(value):
    """Split a comma-separated request value into its non-empty, stripped parts."""
    return [part.strip() for part in (value or "").split(",") if part.strip()]


def _positive_int(value, name):
    if not str(value).isdigit() or int(value) < 1:
        raise ValueError(f"{name} must be a whole number greater than 0")
    return int(value)


def parse_criteria(args):
    """
    Build Criteria from request arguments (a dict of strings). Raises ValueError on bad input.

    Unknown ownership values are kept: the visibility clauses treat a filter with no
    recognized value as matching nothing, the same as the existing search pages.
    """
    types = []
    for kind in _csv(args.get("types")):
        if kind not in TYPES:
            raise ValueError(f"Unknown result type: {kind}")
        if kind not in types:
            types.append(kind)
    if set(types) == set(TYPES):
        types = []      # every type selected is the same as All

    organism_ids = []
    for organism_id in _csv(args.get("organism_ids")):
        if not organism_id.isdigit():
            raise ValueError("Organism IDs must be numbers")
        organism_ids.append(int(organism_id))

    date_added = args.get("date_added") or None
    if date_added is not None and date_added not in DATE_RANGES:
        raise ValueError(f"Date added must be one of: {', '.join(DATE_RANGES)}")

    sort = args.get("sort") or "relevance"
    if sort not in SORTS:
        raise ValueError(f"Sort must be one of: {', '.join(SORTS)}")

    return Criteria(
        q=(args.get("q") or "").strip(),
        types=types,
        organism_ids=organism_ids,
        ownership=[value.lower() for value in _csv(args.get("ownership"))],
        date_added=date_added,
        sort=sort,
        page=_positive_int(args.get("page") or 1, "Page"),
        limit=min(_positive_int(args.get("limit") or DEFAULT_LIMIT, "Limit"), MAX_LIMIT),
    )


def _sort_hits(hits, sort):
    """
    Return hits in display order.

    Python's sort is stable, so the sorts run from the least to the most important key:
    title, then date (newest first, undated last), then the main key for the chosen sort.
    """
    by_title = sorted(hits, key=lambda hit: hit.title.casefold())
    if sort == "title":
        return by_title

    # "" sorts before every date, so reversing puts undated hits last
    newest_first = sorted(by_title, key=lambda hit: hit.date or "", reverse=True)
    if sort == "date":
        return newest_first

    # Relevance. Raw scores from different tables aren't comparable, so each score is scaled
    #  by the best score of its own type. With no keywords every score is 0, which leaves the
    #  newest-first order.
    best = {}
    for hit in hits:
        best[hit.type] = max(best.get(hit.type, 0), hit.score or 0)

    def scaled(hit):
        return (hit.score or 0) / best[hit.type] if best[hit.type] > 0 else 0

    # Hits with the same scaled score take turns by type, so one type with many recent ties
    #  can't fill the first page. A hit's turn is its place among same-type, same-score hits,
    #  counted newest first; within a turn the order stays newest first.
    seen = {}
    turn = {}
    for hit in newest_first:
        group = (hit.type, scaled(hit))
        turn[id(hit)] = seen.get(group, 0)
        seen[group] = turn[id(hit)] + 1

    return sorted(newest_first, key=lambda hit: (-scaled(hit), turn[id(hit)]))


def _pagination(total, page, limit):
    """The same keys the existing search CGIs return."""
    total_pages = ceil(total / limit)
    return {
        "total_results": total,
        "current_page": page,
        "limit": limit,
        "total_pages": total_pages,
        "next_page": page + 1 if page < total_pages else None,
        "prev_page": page - 1 if page > 1 else None,
    }


def _hydrate_page(page_hits, user, cursor, searchers):
    """Hydrate the page with one hydrate() call per type, keeping the page's order."""
    keys_by_type = {}
    for hit in page_hits:
        keys_by_type.setdefault(hit.type, []).append(hit.key)

    hydrated = {kind: searchers[kind].hydrate(keys, user, cursor) for kind, keys in keys_by_type.items()}

    results = []
    for hit in page_hits:
        result = hydrated[hit.type].get(hit.key)
        if result is None:
            # Removed, or made private, between find() and hydrate()
            continue
        result.setdefault("type", hit.type)
        results.append(result)
    return results


def run_search(criteria, user, cursor, searchers):
    """
    Run every searcher and return {"counts", "results", "pagination"}.

    searchers maps each type in TYPES to an object with find() and hydrate().
    Counts cover every type, whatever the type filter, so each filter button can show one.
    """
    hits_by_type = {kind: searchers[kind].find(criteria, user, cursor) for kind in TYPES}

    counts = {kind: len(hits) for kind, hits in hits_by_type.items()}
    counts = {"all": sum(counts.values()), **counts}

    selected = criteria.types or TYPES
    hits = [hit for kind in selected for hit in hits_by_type[kind]]
    hits = _sort_hits(hits, criteria.sort)

    start = (criteria.page - 1) * criteria.limit
    page_hits = hits[start:start + criteria.limit]

    return {
        "counts": counts,
        "results": _hydrate_page(page_hits, user, cursor, searchers),
        "pagination": _pagination(len(hits), criteria.page, criteria.limit),
    }
