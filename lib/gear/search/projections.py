"""
Searcher for saved projections (#1383).

A result is one dataset and one pattern source (a gene list) that was projected onto it. Each
pair can have several saved runs, one per algorithm and z-score setting. The runs are listed in
the result's preview, each linking to the Projection Tool.

Saved runs are read from the registry files the Projection Tool writes:

    <PROJECTIONS_ROOT>/by_dataset/<dataset_id>/projections.json
    {"<gene list share_id>": [{"uuid": "...", "algorithm": "pca", "zscore": false, ...}, ...]}

A run counts only if its result file, <uuid>.csv, is next to the registry. The registry has no
dates, so a pair's date is when its newest result file was written.

The tool behind projections is never named in the UI, but its name still matches in searches.
"""

import json
import sys
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path
from urllib.parse import urlencode

from gear.search import Hit
from gear.search.lookups import visible_datasets, visible_gene_carts
from gear.search.text import compile_terms, score

# TODO #1143: projections are planned to move under www/cache/; change only this path then
PROJECTIONS_ROOT = Path(__file__).resolve().parents[3] / "www" / "projections"

# Older registries key gene lists as "cart.<share_id>"
LEGACY_KEY_PREFIX = "cart."

# How much a search term counts for, by where it was found
TITLE_WEIGHT = 5        # the dataset title or the gene list label
DESCRIPTION_WEIGHT = 3  # their descriptions
CONCEPT_WEIGHT = 1      # words that describe every projection

# Searching for any of these finds every projection. The tool's name is here so it still
#  matches, although it is never shown.
CONCEPT_TEXT = "projection projections projectR pattern patterns"

# The "Date added" choices, as days back from today (as for study pages)
DATE_RANGE_DAYS = {"week": 7, "month": 30, "year": 365}

# Algorithm names as the Projection Tool's menu shortens them
ALGORITHM_LABELS = {"pca": "PCA", "nmf": "NMF", "fixednmf": "Fixed-weight NMF", "binary": "Binary gene count"}


@dataclass(frozen=True)
class SavedRun:
    uuid: str
    algorithm: str
    zscore: bool
    date: str           # YYYY-MM-DD the result file was written


def _read_registry(registry_path):
    """The registry as a dict, or None (with a message) if it is unreadable or not an object."""
    try:
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as err:
        print(f"projections: skipping {registry_path}: {err}", file=sys.stderr)
        return None
    return registry if isinstance(registry, dict) else None


def _saved_run(folder, config):
    """A SavedRun for one registry entry, or None if it is malformed or has no result file."""
    if not isinstance(config, dict):
        return None
    uuid = str(config.get("uuid") or "")
    # The uuid names a file in this folder, so it must not be a path
    if not uuid or Path(uuid).name != uuid:
        return None

    csv_path = folder / f"{uuid}.csv"
    if not csv_path.is_file():
        return None

    # Legacy entries say is_pca instead of algorithm, and have no zscore (read as the Projection Tool does)
    algorithm = config.get("algorithm")
    if "is_pca" in config:
        algorithm = "pca" if config["is_pca"] == 1 else "nmf"
    if not algorithm:
        return None

    written = date.fromtimestamp(csv_path.stat().st_mtime).isoformat()
    return SavedRun(uuid, str(algorithm), bool(config.get("zscore", False)), written)


def read_saved_runs(root=PROJECTIONS_ROOT):
    """
    Every saved run under root, grouped by pair: {(dataset_id, gene list share_id): [SavedRun]}.
    Runs keep their registry order, and each algorithm and z-score setting appears once.
    """
    saved = {}
    for registry_path in sorted(Path(root).glob("by_dataset/*/projections.json")):
        registry = _read_registry(registry_path)
        if not registry:
            continue

        folder = registry_path.parent
        for key, configs in registry.items():
            if not isinstance(configs, list):
                continue
            cart_id = key.removeprefix(LEGACY_KEY_PREFIX)
            runs = saved.setdefault((folder.name, cart_id), [])
            for config in configs:
                run = _saved_run(folder, config)
                if run and all((run.algorithm, run.zscore) != (r.algorithm, r.zscore) for r in runs):
                    runs.append(run)

    return {pair: runs for pair, runs in saved.items() if runs}


def _pair_key(dataset_id, cart_id):
    return f"{dataset_id}/{cart_id}"


def _run_url(dataset, cart_id, run):
    params = {"share_id": dataset["share_id"], "projection_source": cart_id,
              "projection_algorithm": run.algorithm, "zscore": int(run.zscore)}
    return f"./projection.html?{urlencode(params)}"


def _run_label(run):
    label = ALGORITHM_LABELS.get(run.algorithm, run.algorithm)
    return f"{label} · z-scored" if run.zscore else label


class ProjectionSearcher:
    def __init__(self, root=PROJECTIONS_ROOT):
        self.root = Path(root)

    def _visible_pairs(self, runs, user, cursor, ownership=(), organism_ids=()):
        """
        [(dataset, cart, runs)] for the pairs whose dataset and gene list the user may both see.
        The ownership and organism filters narrow the datasets.
        """
        dataset_ids = sorted({dataset_id for dataset_id, _ in runs})
        datasets = {row["id"]: row for row in visible_datasets(cursor, user, ids=dataset_ids,
                                                                ownership=ownership, organism_ids=organism_ids)}
        cart_ids = sorted({cart_id for dataset_id, cart_id in runs if dataset_id in datasets})
        carts = visible_gene_carts(cursor, user, cart_ids) if cart_ids else {}

        return [(datasets[dataset_id], carts[cart_id], pair_runs)
                for (dataset_id, cart_id), pair_runs in runs.items()
                if dataset_id in datasets and cart_id in carts]

    def find(self, criteria, user, cursor):
        runs = read_saved_runs(self.root)

        if criteria.date_added:
            since = (date.today() - timedelta(days=DATE_RANGE_DAYS[criteria.date_added])).isoformat()
            runs = {pair: pair_runs for pair, pair_runs in runs.items()
                    if max(run.date for run in pair_runs) >= since}

        if not runs:
            return []

        patterns = compile_terms(criteria.q)
        hits = []
        for dataset, cart, pair_runs in self._visible_pairs(runs, user, cursor, criteria.ownership,
                                                            criteria.organism_ids):
            pair_score = 0
            if patterns:
                pair_score = score(patterns, [
                    (dataset["title"], TITLE_WEIGHT), (cart["label"], TITLE_WEIGHT),
                    (dataset["ldesc"], DESCRIPTION_WEIGHT), (cart["ldesc"], DESCRIPTION_WEIGHT),
                    (CONCEPT_TEXT, CONCEPT_WEIGHT),
                ])
                if not pair_score:
                    continue
            newest = max(run.date for run in pair_runs)
            hits.append(Hit("projection", _pair_key(dataset["id"], cart["share_id"]), pair_score,
                            dataset["title"], newest))
        return hits

    def hydrate(self, keys, user, cursor):
        runs = read_saved_runs(self.root)
        wanted = {}
        for key in keys:
            dataset_id, _, cart_id = key.partition("/")
            if (dataset_id, cart_id) in runs:
                wanted[(dataset_id, cart_id)] = runs[(dataset_id, cart_id)]
        if not wanted:
            return {}

        results = {}
        for dataset, cart, pair_runs in self._visible_pairs(wanted, user, cursor):
            result = self._result(dataset, cart, pair_runs)
            results[result["id"]] = result
        return results

    def _result(self, dataset, cart, pair_runs):
        cart_id = cart["share_id"]
        runs = [{"label": _run_label(run), "url": _run_url(dataset, cart_id, run)} for run in pair_runs]
        newest = max(run.date for run in pair_runs)
        run_count = f"{len(runs)} saved run{'s' if len(runs) != 1 else ''}"

        meta = [["Pattern source", cart["label"]], ["Organism", dataset["organism"]],
                ["Owner", dataset["owner_name"]], ["Last run", newest]]

        return {
            "type": "projection",
            "id": _pair_key(dataset["id"], cart_id),
            "share_id": dataset["share_id"],
            "title": dataset["title"],
            "desc": dataset["ldesc"],
            "access": dataset["access"],
            "organism": dataset["organism"],
            "date": newest,
            "subline": f"Projection · {cart['label']} · {run_count}",
            "meta": [row for row in meta if row[1]],
            "chips": [dataset["organism"]] if dataset["organism"] else [],
            "links": [],
            "actions": {
                "primary": {"label": "Open in Projection Tool", "url": runs[0]["url"]},
                "secondary": {"label": "Open pattern source in Projection Tool",
                              "url": f"./projection.html?{urlencode({'projection_source': cart_id})}"},
                "share_url": runs[0]["url"],
            },
            "extra": {"runs": runs},
        }
