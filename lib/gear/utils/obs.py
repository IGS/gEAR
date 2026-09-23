# obs.py - obs-dataframe sanitization/categorization helpers for AnnData uploads.

import pandas as pd


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
