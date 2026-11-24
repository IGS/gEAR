"""
Searcher for gene lists (#1383).

Matching follows www/cgi/search_gene_carts.cgi: MySQL full-text search on the label and
description, plus LIKE on the label so words inside underscore names are found. Visibility
follows the same rules (gear.search.access.gene_cart_visibility).

Results preview the first genes of each list. Unweighted lists keep their genes in the
gene_cart_member table; weighted lists keep them in www/carts/cart.<share_id>.tab.
"""

import csv
import sys
from pathlib import Path
from urllib.parse import urlencode

from gear.search import Hit
from gear.search.access import gene_cart_visibility
from gear.search.fulltext import to_boolean_mode_query, to_like_patterns
from gear.search.lookups import access_label

CARTS_DIR = Path(__file__).resolve().parents[3] / "www" / "carts"

# How many genes a result previews
GENE_PREVIEW_COUNT = 12

# The "Date added" choices as MySQL intervals. Only these fixed strings are put in the query.
DATE_INTERVALS = {"week": "INTERVAL 1 WEEK", "month": "INTERVAL 1 MONTH", "year": "INTERVAL 1 YEAR"}

# A list found only by LIKE (e.g. a word inside an underscore name) has no full-text score.
#  This keeps it above zero, below any full-text match.
LIKE_ONLY_SCORE = 0.01


def _placeholders(values):
    return ", ".join(["%s"] * len(values))


def _iso_date(value):
    return value.strftime("%Y-%m-%d") if value else None


def _display_date(value):
    return f"{value:%b} {value.day}, {value.year}" if value else None


class GeneListSearcher:
    def __init__(self, carts_dir=CARTS_DIR):
        self.carts_dir = Path(carts_dir)

    def find(self, criteria, user, cursor):
        visibility, params = gene_cart_visibility(user, criteria.ownership)
        conditions = [visibility]

        fulltext_query = to_boolean_mode_query(criteria.q, prefix=True, boost_phrase=True)
        score_sql = "0"
        score_params = []
        if fulltext_query:
            like_patterns = to_like_patterns(criteria.q)
            score_sql = "MATCH(gc.label, gc.ldesc) AGAINST(%s IN BOOLEAN MODE)"
            score_params = [fulltext_query]
            conditions.append("(" + " OR ".join([score_sql] + ["gc.label LIKE %s"] * len(like_patterns)) + ")")
            params.extend([fulltext_query, *like_patterns])

        if criteria.date_added:
            conditions.append(f"gc.date_added >= DATE_SUB(NOW(), {DATE_INTERVALS[criteria.date_added]})")

        if criteria.organism_ids:
            conditions.append(f"gc.organism_id IN ({_placeholders(criteria.organism_ids)})")
            params.extend(criteria.organism_ids)

        cursor.execute(
            f"SELECT gc.share_id, gc.label, gc.date_added, {score_sql} AS rscore"
            f" FROM gene_cart gc WHERE {' AND '.join(conditions)}",
            [*score_params, *params],
        )

        hits = []
        for share_id, label, date_added, rscore in cursor.fetchall():
            score = float(rscore or 0)
            if fulltext_query and not score:
                score = LIKE_ONLY_SCORE
            hits.append(Hit("genelist", share_id, score, label, _iso_date(date_added)))
        return hits

    def hydrate(self, keys, user, cursor):
        keys = list(dict.fromkeys(keys))
        if not keys:
            return {}

        visibility, params = gene_cart_visibility(user, [])
        cursor.execute(
            f"""
            SELECT gc.id, gc.share_id, gc.label, gc.ldesc, gc.gctype, o.label, gc.user_id, g.user_name,
                   gc.is_public, gc.date_added
              FROM gene_cart gc
                   JOIN organism o ON o.id = gc.organism_id
                   LEFT JOIN guser g ON g.id = gc.user_id
             WHERE gc.share_id IN ({_placeholders(keys)}) AND {visibility}
            """,
            [*keys, *params],
        )
        rows = cursor.fetchall()
        genes = self._unweighted_genes(cursor, [row[0] for row in rows if row[4] != "weighted-list"])

        results = {}
        for (cart_id, share_id, label, ldesc, gctype, organism, owner_id, owner_name, is_public,
             date_added) in rows:
            if gctype == "weighted-list":
                cart_genes = self._weighted_genes(share_id)
            else:
                cart_genes = genes.get(cart_id, [])
            results[share_id] = self._result(share_id, label, ldesc, gctype, organism, owner_name, date_added,
                                             access_label(user, owner_id, is_public, False), cart_genes)
        return results

    def _unweighted_genes(self, cursor, cart_ids):
        """{cart id: [gene symbol]} for unweighted lists, from one query for the whole page."""
        if not cart_ids:
            return {}
        cursor.execute(
            f"SELECT gene_cart_id, gene_symbol FROM gene_cart_member"
            f" WHERE gene_cart_id IN ({_placeholders(cart_ids)}) ORDER BY gene_cart_id, id",
            cart_ids,
        )
        genes = {}
        for cart_id, gene_symbol in cursor.fetchall():
            genes.setdefault(cart_id, []).append(gene_symbol)
        return genes

    def _weighted_genes(self, share_id):
        """Gene symbols of a weighted list (its file's second column), or None if the file is unreadable."""
        path = self.carts_dir / f"cart.{share_id}.tab"
        try:
            with open(path, newline="", encoding="utf-8") as cart_file:
                reader = csv.reader(cart_file, delimiter="\t")
                next(reader, None)      # header
                return [row[1] for row in reader if len(row) > 1]
        except OSError as err:
            print(f"gene_lists: no genes for weighted list {share_id}: {err}", file=sys.stderr)
            return None

    def _result(self, share_id, label, ldesc, gctype, organism, owner_name, date_added, access, genes):
        gene_count = len(genes) if genes is not None else None
        kind = "Weighted gene list" if gctype == "weighted-list" else "Gene list"
        count_text = f"{gene_count} gene{'s' if gene_count != 1 else ''}" if gene_count is not None else None
        url = f"./p?{urlencode({'c': share_id})}"

        meta = [["Genes", str(gene_count) if gene_count is not None else None], ["Organism", organism],
                ["Owner", owner_name], ["Added", _display_date(date_added)]]

        return {
            "type": "genelist",
            "id": share_id,
            "share_id": share_id,
            "title": label,
            "desc": ldesc,
            "access": access,
            "organism": organism,
            "date": _iso_date(date_added),
            "subline": " · ".join(part for part in (kind, organism, count_text) if part),
            "meta": [row for row in meta if row[1]],
            "chips": [organism] if organism else [],
            "links": [],
            "actions": {
                "primary": {"label": "View in Gene Expression", "url": url},
                "secondary": {"label": "Open in Gene List Manager",
                              "url": f"./gene_list_manager.html?{urlencode({'share_id': share_id})}"},
                "share_url": url,
            },
            "extra": {"genes": (genes or [])[:GENE_PREVIEW_COUNT], "gene_count": gene_count},
        }
