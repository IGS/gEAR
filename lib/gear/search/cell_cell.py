"""
Searcher for cell-cell communication analyses (#1383).

Each analysis is a .parquet file in www/cache/domino_signal, with an optional .json file of the
same name describing it (Title, Author, Year, DOI, PMID, "gEAR Dataset URL", Description; read
by gear.dominosignal).

An analysis is linked to a dataset when its "gEAR Dataset URL" names one (p?s=<share_id>,
share_id=<share_id> or dataset_id=<id>). A linked analysis follows that dataset's access: if the
user can't see the dataset, or it no longer exists, the analysis isn't listed. An analysis with
no dataset link counts as public.

The tool that produces these files is never named in the UI, but its name still matches in
searches.
"""

from dataclasses import dataclass
from pathlib import Path
from urllib.parse import parse_qs, quote, urlsplit

from gear import dominosignal
from gear.search import Hit
from gear.search.lookups import visible_datasets
from gear.search.text import compile_terms, score

UNTITLED = "Untitled cell-cell communication analysis"

# How much a search term counts for, by where it was found
TITLE_WEIGHT = 5
DESCRIPTION_WEIGHT = 3
DETAIL_WEIGHT = 1       # author, year, DOI, PMID and the linked dataset's title
CONCEPT_WEIGHT = 1      # words that describe every analysis

# Searching for any of these finds every analysis. The tool's name is here so it still
#  matches, although it is never shown.
CONCEPT_TEXT = ("cell-cell communication dominoSignal ligand ligands receptor receptors "
                "signaling signalling interactions")

# URL parameters that name a dataset, in the order they are checked
_DATASET_PARAMS = (("s", "share_id"), ("share_id", "share_id"), ("dataset_id", "id"))


@dataclass(frozen=True)
class Analysis:
    filename: str           # the .parquet file name
    info: dict | None       # gear.dominosignal's dataset info, or None without a .json file

    def field(self, name):
        """One value from the description, or None."""
        return (self.info or {}).get(name)

    @property
    def title(self):
        return self.field("title") or UNTITLED


def read_analyses(domino_dir=None):
    """Every analysis in domino_dir (default: gear.dominosignal.DOMINO_DIR), by file name."""
    folder = Path(domino_dir or dominosignal.DOMINO_DIR)
    if not folder.is_dir():
        return []
    return [Analysis(path.name, dominosignal.read_dataset_info(str(path.with_suffix(".json"))))
            for path in sorted(folder.glob("*.parquet"))]


def dataset_ref_from_url(url):
    """("share_id", value) or ("id", value) for the dataset a "gEAR Dataset URL" names, else None."""
    if not url:
        return None
    query = parse_qs(urlsplit(url.strip()).query)
    for param, kind in _DATASET_PARAMS:
        values = [value.strip() for value in query.get(param, []) if value.strip()]
        if values:
            return kind, values[0]
    return None


def _linked_datasets(analyses, user, cursor, ownership=(), organism_ids=()):
    """
    {filename: dataset} for the analyses whose linked dataset the user may see, and
    {filename: True} for those with a link at all.
    """
    refs = {}
    for analysis in analyses:
        ref = dataset_ref_from_url(analysis.field("gear_url"))
        if ref:
            refs[analysis.filename] = ref
    if not refs:
        return {}, refs

    ids = sorted({value for kind, value in refs.values() if kind == "id"})
    share_ids = sorted({value for kind, value in refs.values() if kind == "share_id"})
    rows = visible_datasets(cursor, user, ids=ids, share_ids=share_ids, ownership=ownership,
                            organism_ids=organism_ids)
    by_ref = {}
    for row in rows:
        by_ref[("id", row["id"])] = row
        by_ref[("share_id", row["share_id"])] = row

    datasets = {filename: by_ref[ref] for filename, ref in refs.items() if ref in by_ref}
    return datasets, refs


def _doi_url(doi):
    for prefix in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "doi:"):
        if doi.lower().startswith(prefix):
            doi = doi[len(prefix):]
            break
    return f"https://doi.org/{doi.strip()}"


class CellCellSearcher:
    def __init__(self, domino_dir=None):
        self.domino_dir = domino_dir

    def find(self, criteria, user, cursor):
        # The analyses have no date, so a "Date added" filter leaves them all out
        if criteria.date_added:
            return []

        analyses = read_analyses(self.domino_dir)
        datasets, refs = _linked_datasets(analyses, user, cursor, criteria.ownership, criteria.organism_ids)

        # Unlinked analyses count as public and have no organism. Logged-out searches ignore the
        #  ownership filter, as they do for datasets.
        unlinked_shown = (not criteria.organism_ids
                          and (not user or not criteria.ownership or "public" in criteria.ownership))

        patterns = compile_terms(criteria.q)
        hits = []
        for analysis in analyses:
            dataset = datasets.get(analysis.filename)
            if analysis.filename in refs and not dataset:
                continue
            if analysis.filename not in refs and not unlinked_shown:
                continue

            analysis_score = 0
            if patterns:
                details = [analysis.field(name) for name in ("author", "year", "doi", "pmid")]
                if dataset:
                    details.append(dataset["title"])
                analysis_score = score(patterns, [
                    (analysis.field("title"), TITLE_WEIGHT),
                    (analysis.field("description"), DESCRIPTION_WEIGHT),
                    (" ".join(detail for detail in details if detail), DETAIL_WEIGHT),
                    (CONCEPT_TEXT, CONCEPT_WEIGHT),
                ])
                if not analysis_score:
                    continue
            hits.append(Hit("ccc", analysis.filename, analysis_score, analysis.title, None))
        return hits

    def hydrate(self, keys, user, cursor):
        analyses = {analysis.filename: analysis for analysis in read_analyses(self.domino_dir)}
        wanted = [analyses[key] for key in keys if key in analyses]
        if not wanted:
            return {}

        datasets, refs = _linked_datasets(wanted, user, cursor)
        return {analysis.filename: self._result(analysis, datasets.get(analysis.filename))
                for analysis in wanted
                # A linked analysis whose dataset the user can't see is left out
                if analysis.filename not in refs or analysis.filename in datasets}

    def _result(self, analysis, dataset):
        url = f"./cell_cell.html?domino_file={quote(analysis.filename, safe='')}"
        organism = dataset["organism"] if dataset else None
        author = analysis.field("author")
        year = analysis.field("year")
        pmid = analysis.field("pmid") or ""
        doi = analysis.field("doi")

        citation = ", ".join(part for part in (author, year) if part)
        meta = [["Dataset", dataset["title"] if dataset else None], ["Organism", organism],
                ["Authors", author], ["Year", year], ["PMID", pmid]]

        links = []
        if doi:
            links.append({"label": "Read paper", "url": _doi_url(doi)})
        if pmid.isdigit():
            links.append({"label": "PubMed", "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})

        secondary = None
        if dataset:
            secondary = {"label": "View in Gene Expression", "url": f"./p?s={quote(dataset['share_id'], safe='')}"}

        return {
            "type": "ccc",
            "id": analysis.filename,
            "share_id": None,
            "title": analysis.title,
            "desc": analysis.field("description"),
            "access": dataset["access"] if dataset else "Public",
            "organism": organism,
            "date": None,
            "subline": " · ".join(part for part in ("Cell-cell communication", organism, citation) if part),
            "meta": [row for row in meta if row[1]],
            "chips": [organism] if organism else [],
            "links": links,
            "actions": {"primary": {"label": "View signaling results", "url": url},
                        "secondary": secondary,
                        "share_url": url},
            "extra": {"dataset": {"title": dataset["title"], "url": secondary["url"]} if dataset else None},
        }
