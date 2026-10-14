"""Regression tests for closed GitHub issues in gear.metadata (saving upload metadata to MySQL)."""

import json
import re
from unittest import mock

import numpy as np
import pandas as pd
import pytest
from gear import metadata as md

BASE_METADATA = {
    "dataset_uid": "DS1",
    "share_uid": "SHARE1",
    "owner_id": 1,
    "title": "A dataset",
    "sample_taxid": "10090",
    "summary": "Long description",
    "dataset_type": "single-cell",
    "user_pii_affirmed": 1,
}


@pytest.fixture
def cursor(monkeypatch) -> mock.MagicMock:
    """A fake cursor handed out by geardb.Connection in gear.metadata; the taxon always resolves to organism 1."""
    cursor = mock.MagicMock(name="cursor")
    connection = mock.MagicMock(name="connection")
    connection.get_cursor.return_value = cursor
    monkeypatch.setattr(md.geardb, "Connection", lambda: connection)
    monkeypatch.setattr(md.geardb, "get_organism_id_by_taxon_id", lambda taxid: 1)
    return cursor


def save(tmp_path, cursor, **fields) -> dict:
    """Write a metadata.json as store_expression_metadata.cgi does, save it and return the dataset row's values by column."""
    path = tmp_path / "metadata.json"
    path.write_text(json.dumps({**BASE_METADATA, **fields}))
    md.Metadata(file_path=str(path)).save_to_mysql(status="completed")

    (sql, params), = [c.args for c in cursor.execute.call_args_list if "INTO dataset " in c.args[0]]
    columns = [c.strip() for c in re.search(r"dataset \(([^)]*)\)", sql).group(1).split(",")]
    columns.remove("date_added")   # NOW() in the VALUES list, not a parameter
    assert len(columns) == len(params)
    return dict(zip(columns, params))


# ---------------------------------------------------------------- GEO and PubMed IDs


def test_pubmed_without_geo_saves(tmp_path, cursor):
    """Regression for #1113: a PubMed ID with no GEO ID saves (it used to call .strip() on the missing GEO ID)."""
    row = save(tmp_path, cursor, pubmed_id="12345678", geo_accession=None)
    assert row["pubmed_id"] == "12345678"
    assert row["geo_id"] is None


@pytest.mark.parametrize("pubmed_id", ["NA", "n/a", "N/A", "none", "None", "", "  "])
def test_placeholder_pubmed_is_null(tmp_path, cursor, pubmed_id):
    """Regression for #1125: "NA", "n/a" and similar placeholders in the PubMed field save as NULL instead of erroring."""
    row = save(tmp_path, cursor, pubmed_id=pubmed_id)
    assert row["pubmed_id"] is None
    assert row["ldesc"] == "Long description"


@pytest.mark.parametrize("geo_accession", [None, "", "   ", "None", " None "],
                         ids=["null", "empty", "blank", "None-string", "padded-None-string"])
def test_missing_geo_is_null(tmp_path, cursor, geo_accession):
    """Regression for #1106: an empty, blank or "None" GEO ID saves as NULL, so no GEO link is shown."""
    row = save(tmp_path, cursor, geo_accession=geo_accession)
    assert row["geo_id"] is None


def test_geo_id_is_stripped(tmp_path, cursor):
    """Regression for #1106: a real GEO ID is kept, without surrounding whitespace."""
    assert save(tmp_path, cursor, geo_accession=" GSE12345 ")["geo_id"] == "GSE12345"


@pytest.mark.parametrize("pubmed_id", ["11111111, 22222222 33333333", "11111111,22222222,33333333"])
def test_multiple_pubmed_ids_keep_first_and_append_rest(tmp_path, cursor, pubmed_id):
    """Regression for #1106: with several PubMed IDs the first is saved and the rest are appended to the long description."""
    row = save(tmp_path, cursor, pubmed_id=pubmed_id)
    assert row["pubmed_id"] == "11111111"
    assert row["ldesc"] == "Long description<br>Additional Pubmed IDS: 22222222, 33333333"


def test_multiple_pubmed_ids_without_summary(tmp_path, cursor):
    """Regression for #1106: extra PubMed IDs are appended even when there is no long description."""
    row = save(tmp_path, cursor, pubmed_id="11111111 22222222", summary=None)
    assert row["pubmed_id"] == "11111111"
    assert row["ldesc"] == "<br>Additional Pubmed IDS: 22222222"


# ---------------------------------------------------------------- helpers used by save_to_mysql


@pytest.mark.parametrize("value, expected", [
    (np.nan, True), (float("nan"), True), ("nan", True), ("", True),
    ("GSE1", False), (0, False), ("NA", False),
])
def test_is_na(value, expected):
    """Regression for #1106: is_na treats NaN and empty cells as missing (the basis for storing NULL)."""
    assert md.is_na(value) is expected


def test_get_value_from_df():
    """Regression for #1106: get_value_from_df returns None for missing or empty fields and unwraps {"value": ...} cells."""
    df = pd.DataFrame(
        {"value": [np.nan, "", {"value": "GSE1"}, "GSE2", "Single cell RNA-Seq"]},
        index=["geo_accession", "pubmed_id", "wrapped", "plain", "dataset_type"],
    )
    assert md.get_value_from_df(df, "geo_accession") is None
    assert md.get_value_from_df(df, "pubmed_id") is None
    assert md.get_value_from_df(df, "not_there") is None
    assert md.get_value_from_df(df, "wrapped") == "GSE1"
    assert md.get_value_from_df(df, "plain") == "GSE2"
    assert md.get_value_from_df(df, "dataset_type") == "single-cell-rnaseq"
