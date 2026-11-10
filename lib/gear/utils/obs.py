# obs.py - obs-dataframe sanitization/categorization helpers for AnnData uploads.

"""
obs.py - obs dataframe sanitization and categorization helpers for AnnData uploads, and the
obs column names recognized as an upload's own clusters and embeddings.
"""

import re

import pandas as pd


# Names of obs columns holding an upload's own primary analysis (clusters and embeddings), shared
# by gear.primary_analysis and the workbench's get_embedded_tsne_display.cgi.

# Checked in order; the first one found is the dataset's clustering (gear.primary_analysis copies
# it to obs['louvain']). Includes the sample metadata template's cell_clusters, and the names
# scanpy (leiden) and Seurat (seurat_clusters) use.
# Seurat's orig.ident is left out: it holds the sample or project name, not a cell grouping.
VALID_CLUSTER_COLUMN_NAMES = ['cluster', 'cell_type', 'cluster_label', 'subclass_label', 'joint_cluster_round4_annot',
                              'cell_clusters', 'clusters', 'seurat_clusters', 'leiden', 'louvain']
# Seurat's FindClusters also writes one column per resolution, e.g. RNA_snn_res.0.8
SEURAT_CLUSTER_COLUMN_PATTERN = re.compile(r".+_snn_res\.[\d.]+")

# Lowercase pairs are what the Seurat (RDS) uploader writes, from reductions named tsne, umap and pca
VALID_TSNE_PAIRS = [['tSNE_1', 'tSNE_2'], ['tSNE1', 'tSNE2'], ['tsne1_combined', 'tsne2_combined'], ['tsne_1', 'tsne_2'],
                    # These are all carlo's custom ones.  Need to resolve this a different way later
                    ['PC1%var6.14', 'PC2%var1.79']]
VALID_UMAP_PAIRS = [['uMAP_1', 'uMAP_2'], ['uMAP1', 'uMAP2'],
                    ['UMAP_1', 'UMAP_2'], ['UMAP1', 'UMAP2'], ['umap_1', 'umap_2']]
VALID_PCA_PAIRS = [['PCA_1', 'PCA_2'], ['PC_1', 'PC_2'], ['PC1', 'PC2'], ['pca_1', 'pca_2']]


def find_cluster_column(cols: list) -> str | None:
    """
    Return the first column in VALID_CLUSTER_COLUMN_NAMES found in cols, else the first Seurat
    resolution column (e.g. RNA_snn_res.0.8), else None.
    """
    for vname in VALID_CLUSTER_COLUMN_NAMES:
        if vname in cols:
            return vname

    for col in cols:
        if SEURAT_CLUSTER_COLUMN_PATTERN.fullmatch(str(col)):
            return col

    return None


def find_column_pair(cols: list, pairs: list) -> list | None:
    """
    Return the first [x, y] pair from pairs (e.g. VALID_TSNE_PAIRS) with both columns in cols, else None.
    """
    for pair in pairs:
        if pair[0] in cols and pair[1] in cols:
            return pair

    return None


def flag_ambiguous_obs_columns(obs: pd.DataFrame, max_unique: int = 30) -> dict:
    """
    Identify obs columns whose current dtype looks questionable: numeric,
    but with few enough unique values that they're plausibly mislabeled
    categorical data (e.g. replicate numbers, slide numbers stored as
    ints/floats).

    Used by the uploader's "review column types" step so a user can
    explicitly confirm or override how such columns should be treated,
    rather than guessing via a fixed cardinality threshold at display/plot
    time. Note this runs on the obs table *after* categorize_observation_columns()
    has already applied its fixed name-based rules (see anndata_processor.py),
    so it also gives the user a chance to override those hardcoded calls
    (e.g. 'replicate') if they disagree, not just catch columns the fixed
    list misses.

    Parameters
    ----------
    obs : pd.DataFrame
        The obs dataframe to inspect.
    max_unique : int, optional (default: 30)
        Columns with more unique values than this are assumed to be
        genuinely continuous and are not flagged.

    Returns
    -------
    dict
        Keyed by column name, e.g.:
            {
              "replicate": {
                  "current_dtype": "int64",
                  "n_unique": 3,
                  "sample_values": [1, 2, 3],
                  "suggested_type": "categorical",
              }
            }
        Columns that are already categorical/string, or numeric with more
        than `max_unique` unique values, are omitted entirely -- there's
        nothing worth asking the user about for those.
    """

    max_sample_values = 6
    questionable = {}

    for col in obs.columns:
        series = obs[col]

        if isinstance(series.dtype, pd.CategoricalDtype) or series.dtype == object:
            continue
        if not pd.api.types.is_numeric_dtype(series):
            continue

        n_unique = int(series.nunique(dropna=True))
        has_many_uniques = n_unique == 0 or n_unique > max_unique

        # If the column name doesn't contain "cluster" and it has many unique values, we skip it.
        if "cluster" not in col and has_many_uniques:
            continue

        sample_values = sorted(series.dropna().unique().tolist())[:max_sample_values]
        questionable[col] = {
            "current_dtype": str(series.dtype),
            "n_unique": n_unique,
            "sample_values": sample_values,
            # Default guess to pre-select in the review UI -- never applied
            # on its own, the user always makes the final call.
            "suggested_type": "categorical" if n_unique <= 10 else "continuous",
        }

    return questionable


def is_categorical_obs_column(series: pd.Series) -> bool:
    """True if an obs column holds categories (categorical, string or boolean) rather than numbers."""
    return (isinstance(series.dtype, pd.CategoricalDtype) or series.dtype == object
            or pd.api.types.is_bool_dtype(series) or pd.api.types.is_string_dtype(series))


def find_primary_variable_candidates(obs: pd.DataFrame, max_unique: int = 200) -> dict:
    """
    Find the obs columns a user could choose as the dataset's primary variable (#435): the column
    the dataset is mainly about, like a condition or cell type, used to build default displays.

    Only categorical columns with 2 to max_unique values are offered. Columns with one value,
    or a different value for every observation (like cell barcodes), can't group anything.
    Numeric columns are left out; the uploader offers the ambiguous ones from
    flag_ambiguous_obs_columns when the user marks them categorical.

    Returns a dict keyed by column name, e.g.
        {"condition": {"n_unique": 2, "sample_values": ["control", "treated"]}}
    """
    max_sample_values = 6
    candidates = {}

    for col in obs.columns:
        series = obs[col]
        if not is_categorical_obs_column(series):
            continue

        n_unique = int(series.nunique(dropna=True))
        if n_unique < 2 or n_unique > max_unique or n_unique >= len(obs):
            continue

        values = sorted(str(v) for v in series.dropna().unique())
        candidates[str(col)] = {"n_unique": n_unique, "sample_values": values[:max_sample_values]}

    return candidates


def apply_obs_dtype_choices(obs: pd.DataFrame, choices: dict) -> pd.DataFrame:
    """
    Apply user-chosen dtypes to obs columns flagged by flag_ambiguous_obs_columns.

    Parameters
    ----------
    obs : pd.DataFrame
        The obs dataframe to update.
    choices : dict
        Mapping of {column_name: "categorical" | "continuous"}.

    Returns
    -------
    pd.DataFrame
        The same obs dataframe, with the requested columns retyped in
        place. Columns named in `choices` that aren't present in `obs`
        are silently skipped (the file may have changed shape since the
        columns were flagged).
    """

    for col, kind in choices.items():
        if col not in obs.columns:
            continue
        if kind == "categorical":
            obs[col] = obs[col].astype(str).astype("category")
        elif kind == "continuous":
            obs[col] = pd.to_numeric(obs[col], errors="coerce")
        else:
            raise ValueError(f"Unknown dtype choice '{kind}' for column '{col}'")

    return obs

def sanitize_obs_for_h5ad(obs_df: pd.DataFrame) -> pd.DataFrame:
    """Sanitize observation dataframe for downstream storage."""
    for col in obs_df.columns:
        # Convert object columns to categorical, filling NaNs with empty strings first
        if obs_df[col].dtype == 'object':
            obs_df[col] = obs_df[col].fillna('').astype(str)
            obs_df[col] = pd.Categorical(obs_df[col])

        # If all numeric float values are actually integers (all end in .0), convert to int type
        elif pd.api.types.is_float_dtype(obs_df[col]):
            if (obs_df[col].dropna() % 1 == 0).all():
                obs_df[col] = obs_df[col].astype('Int64')  # Use nullable integer type
    return obs_df

def categorize_standard_obs_columns(obs_df: pd.DataFrame) -> pd.DataFrame:
    """Categorize and convert specific observation columns."""
    for str_type in ['cell_type', 'condition', 'replicate', 'time_point', 'time_unit']:
        if str_type in obs_df.columns:
            obs_df[str_type] = pd.Categorical(obs_df[str_type])

    for num_type in ['time_point_order']:
        if num_type in obs_df.columns:
            obs_df[num_type] = pd.to_numeric(obs_df[num_type])

    return obs_df

def standardize_and_sanitize_obs(obs_df: pd.DataFrame) -> pd.DataFrame:
    """Apply standard column-type rules, then general sanitization, to an obs dataframe."""
    obs_df = categorize_standard_obs_columns(obs_df)
    return sanitize_obs_for_h5ad(obs_df)
