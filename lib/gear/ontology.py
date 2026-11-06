"""
ontology.py - Ontology term lookups for dataset and sample metadata (GitHub issue #422).

Sample-level metadata gives tissues as Uberon term IDs and cell types as Cell Ontology (CL)
term IDs. After an upload is converted to H5AD/Zarr, annotate_sample_metadata():
  - adds a column of term names next to each term ID column in adata.obs
    (tissue_ontology -> tissue_ontology_name, cell_ontology -> cell_ontology_name)
  - collects the distinct term IDs, which are saved per dataset (dataset_ontology_term) so
    datasets can be searched by tissue and cell type
  - warns about IDs that couldn't be resolved and recommended columns that are missing

Term names come from the EBI Ontology Lookup Service (OLS4) and are cached in the
ontology_term table. Lookups never fail an upload: unresolved terms only produce warnings.
"""

import re
import sys

import geardb
import pandas as pd
import requests
from gear.metadata_fields import SAMPLE_FIELDS, parse_yes_no
from gear.primary_analysis import detect_primary_analysis
from gear.serverconfig import ServerConfig
from gear.utils.obs import (
    VALID_PCA_PAIRS,
    VALID_TSNE_PAIRS,
    VALID_UMAP_PAIRS,
    find_cluster_column,
    find_column_pair,
)

OLS_API_URL = "https://www.ebi.ac.uk/ols4/api"
OLS_TIMEOUT = 10   # seconds

# Standard obs column -> (dataset_ontology_term.field, name column added to obs, expected prefix)
OBS_ONTOLOGY_COLUMNS = {
    "tissue_ontology": ("tissue_type", "tissue_ontology_name", "UBERON"),
    "cell_ontology": ("cell_type", "cell_ontology_name", "CL"),
}

# Other names people use for those columns, e.g. the CELLxGENE schema's *_ontology_term_id
OBS_COLUMN_ALIASES = {
    "tissue_ontology": ("tissue_ontology_id", "tissue_ontology_term_id", "uberon_id", "uberon"),
    "cell_ontology": ("cell_ontology_id", "cell_type_ontology", "cell_type_ontology_id",
                      "cell_type_ontology_term_id", "cl_id"),
}

DEFAULT_RECOMMENDED_OBS_COLUMNS = ("tissue_ontology", "treatment", "age_value", "age_unit", "sex_assigned_at_birth")

TERM_ID_PATTERN = re.compile(r"([A-Za-z]+)[:_](\d+)")

# adata.uns key holding which obs column holds each standard sample metadata field (see find_standard_obs_columns)
OBS_COLUMN_MAP_KEY = "gear_obs_columns"

# Embedding name in the column map -> obs column pairs recognized for it
EMBEDDING_COLUMN_PAIRS = {"umap": VALID_UMAP_PAIRS, "tsne": VALID_TSNE_PAIRS, "pca": VALID_PCA_PAIRS}


def normalize_term_id(value) -> str | None:
    """'uberon_0002227 ' -> 'UBERON:0002227'. Returns None if the value isn't a term ID."""
    m = TERM_ID_PATTERN.fullmatch(str(value).strip())
    if not m:
        return None
    return f"{m.group(1).upper()}:{m.group(2)}"


def split_term_ids(value) -> list[str]:
    """A cell may hold several comma-separated term IDs. Blank cells give an empty list."""
    if value is None or (isinstance(value, float) and value != value):
        return []
    return [part.strip() for part in str(value).split(",") if part.strip() and part.strip().lower() != "nan"]


def get_recommended_obs_columns() -> list[str]:
    """Sample metadata columns this portal asks for, from gear.ini [dataset_uploader] recommended_obs_columns."""
    config = ServerConfig().parse()
    value = config.get("dataset_uploader", "recommended_obs_columns", fallback=None)
    if value is None:
        return list(DEFAULT_RECOMMENDED_OBS_COLUMNS)
    return [c.strip() for c in value.split(",") if c.strip()]


def lookup_terms(term_ids) -> dict:
    """
    Return {term_id: name} for the given term IDs (normalized, e.g. 'CL:0000601'). Names come from
    the ontology_term cache, then OLS. Terms that can't be found map to None.
    """
    term_ids = sorted({t for t in (normalize_term_id(t) for t in term_ids) if t})
    labels = {t: None for t in term_ids}
    if not term_ids:
        return labels

    labels.update(_read_cache(term_ids))

    missing = [t for t in term_ids if labels[t] is None]
    fetched = _fetch_from_ols(missing)
    labels.update(fetched)
    _write_cache(fetched)

    return labels


def _read_cache(term_ids: list[str]) -> dict:
    try:
        cnx = geardb.Connection()
        cursor = cnx.get_cursor()
        placeholders = ", ".join(["%s"] * len(term_ids))
        cursor.execute(f"SELECT term_id, label FROM ontology_term WHERE term_id IN ({placeholders})", tuple(term_ids))
        rows = cursor.fetchall() or []
        cursor.close()
        cnx.close()
    except Exception as e:
        # The cache only saves lookups, so carry on without it
        print(f"Could not read the ontology_term cache: {e}", file=sys.stderr)
        return {}
    return {term_id: label for term_id, label in rows if label}


def _write_cache(labels: dict) -> None:
    found = [(t, t.split(":")[0].lower(), label) for t, label in labels.items() if label]
    if not found:
        return
    try:
        cnx = geardb.Connection()
        cursor = cnx.get_cursor()
        cursor.executemany("""
            INSERT INTO ontology_term (term_id, ontology, label) VALUES (%s, %s, %s)
            ON DUPLICATE KEY UPDATE label = VALUES(label), date_fetched = NOW()
        """, found)
        cnx.commit()
        cursor.close()
        cnx.close()
    except Exception as e:
        print(f"Could not write the ontology_term cache: {e}", file=sys.stderr)


def _fetch_from_ols(term_ids: list[str]) -> dict:
    """Look up term names in OLS4. Stops at the first connection problem rather than waiting on each term."""
    labels = {}
    with requests.Session() as session:
        for term_id in term_ids:
            ontology = term_id.split(":")[0].lower()
            try:
                response = session.get(f"{OLS_API_URL}/ontologies/{ontology}/terms",
                                       params={"obo_id": term_id}, timeout=OLS_TIMEOUT)
            except requests.RequestException as e:
                print(f"OLS lookup failed, skipping remaining terms: {e}", file=sys.stderr)
                break

            # 404 means OLS doesn't know the ontology or term
            if response.status_code != 200:
                continue
            terms = response.json().get("_embedded", {}).get("terms", [])
            if terms and terms[0].get("label"):
                labels[term_id] = terms[0]["label"]
    return labels


def find_obs_column(obs: pd.DataFrame, standard_name: str) -> str | None:
    """Return the obs column holding a standard sample metadata field, matching names case-insensitively."""
    by_lower = {str(c).lower(): c for c in obs.columns}
    for name in (standard_name, *OBS_COLUMN_ALIASES.get(standard_name, ())):
        if name.lower() in by_lower:
            return by_lower[name.lower()]
    return None


def find_standard_obs_columns(obs: pd.DataFrame) -> dict:
    """
    Return which obs column holds each standard sample metadata field that is present, without renaming
    anything, e.g.
        {"tissue_ontology": "Tissue_Ontology_ID", "tissue_ontology_name": "tissue_ontology_name",
         "cell_clusters": "seurat_clusters", "umap": ["UMAP_1", "UMAP_2"]}

    Keys are the sample metadata template's column names (gear.metadata_fields.SAMPLE_FIELDS), the term
    name columns gEAR adds, "cell_clusters" for the cluster column, and "umap", "tsne" and "pca" for
    coordinate column pairs. Saved in adata.uns so later steps (and a future column picker in the
    uploader) can find these columns without guessing again.
    """
    columns = {}
    embedding_fields = {"UMAP_1", "UMAP_2", "cell_clusters"}
    names = [f.name for f in SAMPLE_FIELDS if f.name not in embedding_fields]
    names += [name_column for _, name_column, _ in OBS_ONTOLOGY_COLUMNS.values()]
    for name in names:
        column = find_obs_column(obs, name)
        if column is not None:
            columns[name] = str(column)

    obs_columns = obs.columns.tolist()
    cluster_column = find_cluster_column(obs_columns)
    if cluster_column is not None:
        columns["cell_clusters"] = cluster_column
    for embedding, pairs in EMBEDDING_COLUMN_PAIRS.items():
        pair = find_column_pair(obs_columns, pairs)
        if pair is not None:
            columns[embedding] = list(pair)

    return columns


def annotate_sample_metadata(obs: pd.DataFrame, recommended_columns=None) -> tuple[pd.DataFrame, dict]:
    """
    Add ontology term name columns to an obs table and check it against the sample metadata template.

    Returns (obs, summary) where summary is:
        {
            "ontology_terms": {"tissue_type": [term IDs], "cell_type": [term IDs]},
            "warnings": [messages for the user],
        }
    Only fields with a term ID column in obs appear in ontology_terms. Column names in the warnings are
    wrapped in backticks, which the uploader shows in bold.
    """
    if recommended_columns is None:
        recommended_columns = get_recommended_obs_columns()

    ontology_terms = {}
    warnings = []

    for standard_name, (field, name_column, prefix) in OBS_ONTOLOGY_COLUMNS.items():
        column = find_obs_column(obs, standard_name)
        if column is None:
            continue

        raw_values = obs[column].astype(object).tolist()
        distinct = sorted({part for value in set(raw_values) for part in split_term_ids(value)})

        valid, invalid, wrong_prefix = [], [], []
        for value in distinct:
            term_id = normalize_term_id(value)
            if term_id is None:
                invalid.append(value)
            elif not term_id.startswith(prefix + ":"):
                wrong_prefix.append(value)
            else:
                valid.append(term_id)

        labels = lookup_terms(valid)
        unresolved = [t for t in valid if labels.get(t) is None]
        ontology_terms[field] = sorted(set(valid))

        def term_names(value):
            names = []
            for part in split_term_ids(value):
                term_id = normalize_term_id(part)
                names.append(labels.get(term_id) or part)
            return ", ".join(names)

        names = [term_names(v) for v in raw_values]
        obs[name_column] = pd.Categorical(names)

        if invalid:
            warnings.append(f"Column `{column}`: these values are not {prefix} term IDs (like {prefix}:0000000): {_preview(invalid)}")
        if wrong_prefix:
            warnings.append(f"Column `{column}`: expected {prefix} terms, but found: {_preview(wrong_prefix)}")
        if unresolved:
            warnings.append(f"Column `{column}`: these terms could not be found in the ontology, so their names "
                            f"weren't added: {_preview(unresolved)}")

    lowered = {str(c).lower() for c in obs.columns}
    missing = [c for c in recommended_columns
               if c.lower() not in lowered and not (c in OBS_ONTOLOGY_COLUMNS and find_obs_column(obs, c))]
    if missing:
        warnings.append("These recommended sample metadata columns are missing: {0}.".format(
            ", ".join(f"`{c}`" for c in missing)))

    return obs, {"ontology_terms": ontology_terms, "warnings": warnings}


def _preview(values: list, limit: int = 5) -> str:
    shown = ", ".join(str(v) for v in values[:limit])
    if len(values) > limit:
        shown += f" (and {len(values) - limit} more)"
    return shown


def annotate_upload_obs(adata, metadata: dict) -> None:
    """
    Run annotate_sample_metadata() on the obs table of an upload being processed (adata.obs is
    replaced), recording the results in its metadata.json contents (metadata, updated in place):
        ontology_terms              saved to dataset_ontology_term when the upload is finalized
        sample_metadata_warnings    shown to the user by the uploader
    Also warns when the metadata says primary analysis is included but gEAR's primary analysis
    won't find any clusters or coordinates, and saves find_standard_obs_columns() in
    adata.uns["gear_obs_columns"]. Never raises: a problem here shouldn't fail the upload.
    """
    try:
        adata.obs, summary = annotate_sample_metadata(adata.obs)
    except Exception as e:
        print(f"Could not annotate sample metadata: {e}", file=sys.stderr)
        summary = {"ontology_terms": {}, "warnings": ["Tissue and cell type ontology terms could not be processed. "
                                                      "Please contact the gEAR team if you need them."]}

    try:
        if parse_yes_no(metadata.get("primary_analysis_completed")) and not detect_primary_analysis(adata):
            summary["warnings"].append(
                "The metadata says primary analysis is included, but no cluster column (like `cell_clusters` or "
                "`seurat_clusters`) or UMAP, tSNE or PCA coordinate columns (like `UMAP_1` and `UMAP_2`) were found.")
    except Exception as e:
        print(f"Could not check for primary analysis results: {e}", file=sys.stderr)

    try:
        adata.uns[OBS_COLUMN_MAP_KEY] = find_standard_obs_columns(adata.obs)
    except Exception as e:
        print(f"Could not record the standard obs columns: {e}", file=sys.stderr)

    metadata["ontology_terms"] = summary["ontology_terms"]
    metadata["sample_metadata_warnings"] = summary["warnings"]
