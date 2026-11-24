"""
Searcher for collections (#1383): the "layout" table.

Collections have only a label: no description, date or full-text index. A collection matches
when its label contains a search word, or when one of its datasets matches the dataset
full-text search. Only datasets the user may see count, so a private dataset's text never
reveals a collection.

Visibility follows www/cgi/get_user_layouts.cgi (gear.search.access.layout_visibility).
Collections have no date, so a "Date added" filter leaves them all out.

(The module is named after the table: a module named "collections" would hide Python's own.)
"""

import re
from urllib.parse import urlencode

from gear.search import Hit
from gear.search.access import (
    DATASET_LIVE_CLAUSE,
    dataset_visibility,
    layout_visibility,
)
from gear.search.fulltext import (
    MIN_PARTIAL_MATCH_LENGTH,
    _split_on_quotes,
    _words,
    to_boolean_mode_query,
    to_like_patterns,
)
from gear.search.lookups import access_label
from gear.search.text import compile_terms, score

# How many of a collection's datasets a result lists
DATASET_PREVIEW_COUNT = 10

# A search term found in the label counts as much as a title match elsewhere. Each matching
#  dataset adds 1, up to MAX_DATASET_SCORE, so a label match always ranks first.
LABEL_WEIGHT = 5
MAX_DATASET_SCORE = LABEL_WEIGHT - 1

# The same columns search_datasets.cgi searches
DATASET_MATCH = "MATCH(d.title, d.ldesc, d.geo_id, d.pubmed_id) AGAINST(%s IN BOOLEAN MODE)"

# The member datasets the user may see, joined to each collection
_MEMBER_JOINS = """
       LEFT JOIN layout_displays ld ON ld.layout_id = l.id
       LEFT JOIN dataset_display dd ON dd.id = ld.display_id
       LEFT JOIN dataset d ON d.id = dd.dataset_id AND {live} AND {visible}
"""


def _placeholders(values):
    return ", ".join(["%s"] * len(values))


def _label_prefilter_patterns(text):
    """
    LIKE patterns that let SQL skip collections whose label can't match: one "%word%" per word of
    3 or more characters, quoted or not. Python then checks each label exactly (_label_score).
    """
    words = [word for _, chunk in _split_on_quotes(text) for word in _words(chunk)]
    return list(dict.fromkeys(f"%{word}%" for word in words if len(word) >= MIN_PARTIAL_MATCH_LENGTH))


def _like_regex(pattern):
    """A Python regex for a LIKE pattern from to_like_patterns (case-insensitive, as MySQL's default)."""
    parts = re.split(r"(\\.|%|_)", pattern)
    converted = []
    for part in parts:
        if part == "%":
            converted.append(".*")
        elif part == "_":
            converted.append(".")
        elif part.startswith("\\") and len(part) == 2:
            converted.append(re.escape(part[1]))
        else:
            converted.append(re.escape(part))
    return re.compile("".join(converted), re.IGNORECASE | re.DOTALL)


def _label_score(label, q):
    """
    How well a label matches the search text: the keyword rules every searcher uses, or the
    search pages' LIKE rule for words inside underscore names ("potential" in "Groves_potential").
    """
    label_score = score(compile_terms(q), [(label, LABEL_WEIGHT)])
    if not label_score and any(_like_regex(pattern).fullmatch(label or "") for pattern in to_like_patterns(q)):
        label_score = LABEL_WEIGHT
    return label_score


class CollectionSearcher:
    def find(self, criteria, user, cursor):
        if criteria.date_added:
            return []

        fulltext_query = to_boolean_mode_query(criteria.q, prefix=True, boost_phrase=True)
        dataset_like = to_like_patterns(criteria.q)

        select_params = []
        matching_sql = "0"
        if fulltext_query:
            member_match = "(" + " OR ".join([DATASET_MATCH] + ["d.title LIKE %s"] * len(dataset_like)) + ")"
            matching_sql = f"COUNT(DISTINCT CASE WHEN {member_match} THEN d.id END)"
            select_params = [fulltext_query, *dataset_like]

        dataset_clause, dataset_params = dataset_visibility(user, [])
        layout_clause, layout_params = layout_visibility(user, criteria.ownership)

        having = []
        having_params = []
        if fulltext_query:
            label_patterns = _label_prefilter_patterns(criteria.q)
            having.append("(" + " OR ".join(["matching_datasets > 0"] + ["l.label LIKE %s"] * len(label_patterns))
                          + ")")
            having_params.extend(label_patterns)
        if criteria.organism_ids:
            having.append(f"COUNT(DISTINCT CASE WHEN d.organism_id IN ({_placeholders(criteria.organism_ids)})"
                          " THEN d.id END) > 0")
            having_params.extend(criteria.organism_ids)

        joins = _MEMBER_JOINS.format(live=DATASET_LIVE_CLAUSE, visible=dataset_clause)
        cursor.execute(
            f"""
            SELECT l.share_id, l.label, {matching_sql} AS matching_datasets
              FROM layout l {joins}
             WHERE {layout_clause}
             GROUP BY l.id, l.share_id, l.label
             {"HAVING " + " AND ".join(having) if having else ""}
            """,
            [*select_params, *dataset_params, *layout_params, *having_params],
        )

        hits = []
        for share_id, label, matching_datasets in cursor.fetchall():
            collection_score = 0
            if fulltext_query:
                collection_score = _label_score(label, criteria.q) + min(int(matching_datasets or 0),
                                                                         MAX_DATASET_SCORE)
                if not collection_score:
                    # Passed SQL's rough label check but not the exact one
                    continue
            hits.append(Hit("collection", share_id, collection_score, label or "", None))
        return hits

    def hydrate(self, keys, user, cursor):
        keys = list(dict.fromkeys(keys))
        if not keys:
            return {}

        layout_clause, layout_params = layout_visibility(user, [])
        cursor.execute(
            f"""
            SELECT l.id, l.share_id, l.label, l.user_id, g.user_name, l.is_public, l.is_domain
              FROM layout l
                   LEFT JOIN guser g ON g.id = l.user_id
             WHERE l.share_id IN ({_placeholders(keys)}) AND {layout_clause}
            """,
            [*keys, *layout_params],
        )
        layouts = cursor.fetchall()
        if not layouts:
            return {}

        members = self._members(cursor, user, [row[1] for row in layouts])
        results = {}
        for _, share_id, label, owner_id, owner_name, is_public, is_domain in layouts:
            # Site-highlighted collections are shown to everyone, so they count as public
            access = access_label(user, owner_id, is_public or is_domain, False)
            results[share_id] = self._result(share_id, label, owner_name, access, members.get(share_id, []))
        return results

    def _members(self, cursor, user, share_ids):
        """{collection share_id: [(dataset share_id, title, organism)]}, each dataset once, in grid order."""
        dataset_clause, dataset_params = dataset_visibility(user, [])
        cursor.execute(
            f"""
            SELECT l.share_id, d.share_id, d.title, o.label
              FROM layout l
                   JOIN layout_displays ld ON ld.layout_id = l.id
                   JOIN dataset_display dd ON dd.id = ld.display_id
                   JOIN dataset d ON d.id = dd.dataset_id
                   JOIN organism o ON o.id = d.organism_id
             WHERE l.share_id IN ({_placeholders(share_ids)}) AND {DATASET_LIVE_CLAUSE} AND {dataset_clause}
             ORDER BY l.share_id, ld.grid_position
            """,
            [*share_ids, *dataset_params],
        )
        members = {}
        for layout_share_id, dataset_share_id, title, organism in cursor.fetchall():
            datasets = members.setdefault(layout_share_id, [])
            if all(dataset_share_id != existing[0] for existing in datasets):
                datasets.append((dataset_share_id, title, organism))
        return members

    def _result(self, share_id, label, owner_name, access, datasets):
        organisms = sorted({organism for _, _, organism in datasets if organism})
        count_text = f"{len(datasets)} dataset{'s' if len(datasets) != 1 else ''}"
        url = f"./p?{urlencode({'l': share_id})}"

        meta = [["Datasets", str(len(datasets))], ["Organisms", ", ".join(organisms)], ["Owner", owner_name]]

        return {
            "type": "collection",
            "id": share_id,
            "share_id": share_id,
            "title": label,
            "desc": None,
            "access": access,
            "organism": ", ".join(organisms) or None,
            "date": None,
            "subline": " · ".join(part for part in ("Collection", count_text, ", ".join(organisms)) if part),
            "meta": [row for row in meta if row[1]],
            "chips": organisms,
            "links": [],
            "actions": {"primary": {"label": "View in Gene Expression", "url": url},
                        "secondary": None,
                        "share_url": url},
            "extra": {
                "datasets": [{"title": title, "url": f"./p?{urlencode({'s': dataset_share_id})}"}
                             for dataset_share_id, title, _ in datasets[:DATASET_PREVIEW_COUNT]],
                "dataset_count": len(datasets),
            },
        }
