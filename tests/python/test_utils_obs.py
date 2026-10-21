"""gear.utils.obs: flagging ambiguous numeric obs columns, applying the user's dtype choices, and h5ad-safe sanitizing."""

import numpy as np
import pandas as pd
import pytest
from gear.utils.obs import (
    apply_obs_dtype_choices,
    flag_ambiguous_obs_columns,
    sanitize_obs_for_h5ad,
)


@pytest.fixture
def obs():
    n = 40
    return pd.DataFrame({
        "replicate": [1, 2, 3, 4] * (n // 4),
        "slide": [float(i % 12) for i in range(n)],
        "n_counts": np.arange(n, dtype=float) * 1.5,
        "leiden_cluster": list(range(n)),
        "cell_type": ["a", "b"] * (n // 2),
        "condition": pd.Categorical([1, 2] * (n // 2)),
    })


def test_flags_low_cardinality_numeric_columns(obs):
    flagged = flag_ambiguous_obs_columns(obs)
    assert flagged["replicate"] == {"current_dtype": "int64", "n_unique": 4, "sample_values": [1, 2, 3, 4],
                                    "suggested_type": "categorical"}
    assert flagged["slide"]["n_unique"] == 12
    assert flagged["slide"]["suggested_type"] == "continuous"
    assert flagged["slide"]["sample_values"] == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]


def test_skips_string_categorical_and_high_cardinality_columns(obs):
    flagged = flag_ambiguous_obs_columns(obs)
    assert "cell_type" not in flagged and "condition" not in flagged and "n_counts" not in flagged


def test_cluster_columns_are_flagged_regardless_of_cardinality(obs):
    flagged = flag_ambiguous_obs_columns(obs)
    assert flagged["leiden_cluster"]["n_unique"] == 40
    assert flagged["leiden_cluster"]["suggested_type"] == "continuous"


def test_max_unique_threshold(obs):
    assert "slide" not in flag_ambiguous_obs_columns(obs, max_unique=10)
    assert "n_counts" in flag_ambiguous_obs_columns(obs, max_unique=40)


def test_all_missing_numeric_column_is_not_flagged():
    assert flag_ambiguous_obs_columns(pd.DataFrame({"score": [np.nan, np.nan]})) == {}


def test_apply_dtype_choices(obs):
    result = apply_obs_dtype_choices(obs, {"replicate": "categorical", "slide": "continuous", "gone": "categorical"})
    assert result is obs
    assert isinstance(obs["replicate"].dtype, pd.CategoricalDtype)
    assert sorted(obs["replicate"].cat.categories) == ["1", "2", "3", "4"]
    assert pd.api.types.is_float_dtype(obs["slide"])


def test_apply_continuous_choice_coerces_bad_values_to_nan():
    obs = pd.DataFrame({"dose": ["1", "2.5", "high"]})
    apply_obs_dtype_choices(obs, {"dose": "continuous"})
    assert obs["dose"].tolist()[:2] == [1.0, 2.5] and np.isnan(obs["dose"].iloc[2])


def test_apply_unknown_choice_raises(obs):
    with pytest.raises(ValueError, match="Unknown dtype choice 'ordinal'"):
        apply_obs_dtype_choices(obs, {"replicate": "ordinal"})


@pytest.fixture
def messy_obs():
    return pd.DataFrame({
        "mixed": [1, "a", None],
        "whole": [1.0, 2.0, np.nan],
        "fraction": [1.5, 2.0, 3.0],
        "count": [1, 2, 3],
    }, index=["c1", "c2", "c3"])


def test_sanitize_obs_dtypes(messy_obs):
    obs = sanitize_obs_for_h5ad(messy_obs)
    assert obs["mixed"].tolist() == ["1", "a", ""]
    assert isinstance(obs["mixed"].dtype, pd.CategoricalDtype)
    assert str(obs["whole"].dtype) == "Int64" and obs["whole"].isna().tolist() == [False, False, True]
    assert obs["fraction"].dtype == float
    assert obs["count"].dtype == "int64"


def test_sanitized_obs_round_trips_through_h5ad(tmp_path, messy_obs):
    ad = pytest.importorskip("anndata")
    path = tmp_path / "obs.h5ad"
    with pytest.raises(TypeError):
        ad.AnnData(X=np.zeros((3, 1)), obs=messy_obs.copy()).write_h5ad(path)

    ad.AnnData(X=np.zeros((3, 1)), obs=sanitize_obs_for_h5ad(messy_obs.copy())).write_h5ad(path)
    obs = ad.read_h5ad(path).obs
    assert obs.index.tolist() == ["c1", "c2", "c3"]
    assert obs["mixed"].tolist() == ["1", "a", ""]
    assert obs["whole"].tolist()[:2] == [1, 2] and pd.isna(obs["whole"].iloc[2])
    assert obs["fraction"].tolist() == [1.5, 2.0, 3.0]
