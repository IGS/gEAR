"""
Searcher for datasets (#1383).

Matching follows www/cgi/search_datasets.cgi (the Dataset Explorer): MySQL full-text search on
the title, description, GEO ID and PubMed ID, plus LIKE on the title so words inside underscore
names are found. Visibility follows the same rules (gear.search.access.dataset_visibility), and
only fully loaded datasets that aren't marked for removal are searched.

The preview image is chosen the same way as in the Dataset Explorer.
"""

from datetime import date
from pathlib import Path
from urllib.parse import quote, urlencode

from gear.search import Hit
from gear.search.access import DATASET_LIVE_CLAUSE, dataset_visibility
from gear.search.fulltext import to_boolean_mode_query, to_like_patterns
from gear.search.lookups import visible_datasets

IMAGE_ROOT = Path(__file__).resolve().parents[3] / "www" / "img" / "dataset_previews"
WEB_IMAGE_ROOT = "./img/dataset_previews"

# The "Date added" choices as MySQL intervals. Only these fixed strings are put in the query.
DATE_INTERVALS = {"week": "INTERVAL 1 WEEK", "month": "INTERVAL 1 MONTH", "year": "INTERVAL 1 YEAR"}

# A dataset found only by LIKE (e.g. a word inside an underscore name) has no full-text score.
#  This keeps it above zero, below any full-text match.
LIKE_ONLY_SCORE = 0.01

# Dataset types as the Dataset Explorer's type filter names them. Other types are shown as stored.
DTYPE_LABELS = {
    "bulk-rnaseq": "Bulk RNA-seq",
    "microarray": "Microarray",
    "single-cell-rnaseq": "Single-cell RNA-seq",
    "scRNA-seq": "Single-cell RNA-seq",
    "single-nucleus-rnaseq": "Single-nucleus RNA-seq",
    "gosling": "Epigenetic",
    "spatial": "Spatial transcriptomics",
    "spatial-h5ad": "Spatial transcriptomics",
}

FULLTEXT_MATCH = "MATCH(d.title, d.ldesc, d.geo_id, d.pubmed_id) AGAINST(%s IN BOOLEAN MODE)"


def _placeholders(values):
    return ", ".join(["%s"] * len(values))


def _display_date(iso_date):
    """"2026-01-02" as "Jan 2, 2026"."""
    if not iso_date:
        return None
    value = date.fromisoformat(iso_date)
    return f"{value:%b} {value.day}, {value.year}"


class DatasetSearcher:
    def __init__(self, image_root=IMAGE_ROOT):
        self.image_root = Path(image_root)

    def find(self, criteria, user, cursor):
        visibility, params = dataset_visibility(user, criteria.ownership)
        conditions = [DATASET_LIVE_CLAUSE, visibility]

        fulltext_query = to_boolean_mode_query(criteria.q, prefix=True, boost_phrase=True)
        score_sql = "0"
        score_params = []
        if fulltext_query:
            like_patterns = to_like_patterns(criteria.q)
            score_sql = FULLTEXT_MATCH
            score_params = [fulltext_query]
            conditions.append("(" + " OR ".join([FULLTEXT_MATCH] + ["d.title LIKE %s"] * len(like_patterns)) + ")")
            params.extend([fulltext_query, *like_patterns])

        if criteria.date_added:
            conditions.append(f"d.date_added >= DATE_SUB(NOW(), {DATE_INTERVALS[criteria.date_added]})")

        if criteria.organism_ids:
            conditions.append(f"d.organism_id IN ({_placeholders(criteria.organism_ids)})")
            params.extend(criteria.organism_ids)

        cursor.execute(
            f"SELECT d.share_id, d.title, d.date_added, {score_sql} AS rscore"
            f" FROM dataset d WHERE {' AND '.join(conditions)}",
            [*score_params, *params],
        )

        hits = []
        for share_id, title, date_added, rscore in cursor.fetchall():
            score = float(rscore or 0)
            if fulltext_query and not score:
                score = LIKE_ONLY_SCORE
            iso_date = date_added.strftime("%Y-%m-%d") if date_added else None
            hits.append(Hit("dataset", share_id, score, title, iso_date))
        return hits

    def hydrate(self, keys, user, cursor):
        keys = list(dict.fromkeys(keys))
        if not keys:
            return {}
        return {row["share_id"]: self._result(row) for row in visible_datasets(cursor, user, share_ids=keys)}

    def _preview_image_url(self, dataset_id, dtype):
        """The Dataset Explorer's choice of image, or None when there is none to show."""
        for name in (f"{dataset_id}.default.png", f"{dataset_id}.single.default.png"):
            if (self.image_root / name).exists():
                return f"{WEB_IMAGE_ROOT}/{name}"
        if dtype == "gosling":
            return f"{WEB_IMAGE_ROOT}/gosling.png"
        return None

    def _result(self, row):
        share_id = row["share_id"]
        dtype = DTYPE_LABELS.get(row["dtype"], row["dtype"]) if row["dtype"] else None
        url = f"./p?{urlencode({'s': share_id})}"

        meta = [["Type", dtype], ["Organism", row["organism"]], ["Owner", row["owner_name"]],
                ["Added", _display_date(row["date_added"])], ["GEO ID", row["geo_id"]],
                ["PubMed ID", row["pubmed_id"]]]

        links = []
        if row["geo_id"]:
            links.append({"label": f"GEO: {row['geo_id']}",
                          "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?{urlencode({'acc': row['geo_id']})}"})
        if row["pubmed_id"]:
            links.append({"label": f"PubMed: {row['pubmed_id']}",
                          "url": f"https://pubmed.ncbi.nlm.nih.gov/{quote(str(row['pubmed_id']), safe='')}/"})

        return {
            "type": "dataset",
            "id": share_id,
            "share_id": share_id,
            "title": row["title"],
            "desc": row["ldesc"],
            "access": row["access"],
            "organism": row["organism"],
            "date": row["date_added"],
            "subline": " · ".join(part for part in ("Dataset", dtype, row["organism"]) if part),
            "meta": [item for item in meta if item[1]],
            "chips": [row["organism"]] if row["organism"] else [],
            "links": links,
            "actions": {
                "primary": {"label": "View in Gene Expression", "url": url},
                "secondary": {"label": "Open in Dataset Explorer",
                              "url": f"./dataset_explorer.html?{urlencode({'share_id': share_id})}"},
                "share_url": url,
            },
            "extra": {"preview_image_url": self._preview_image_url(row["id"], row["dtype"])},
        }
