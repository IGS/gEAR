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


# ---------------------------------------------------------------- v2 metadata template (#390)


def test_v2_fields_are_saved(tmp_path, cursor):
    """Prefixed template values are split/stripped into the dataset columns."""
    row = save(tmp_path, cursor, sample_taxid="NCBI:txid10090", reference_annot_id="RefSeq:GCF_000001635.27",
               repository_accession="GEO:GSE168901", pubmed_id="PMID:25855195", device="OBI:0002630",
               assay="OBI:0002631", log_transformation="log2 transformation", primary_analysis_completed="Yes",
               dataset_format="h5ad")
    assert row["annotation_source"] == "RefSeq"
    assert row["annotation_release"] == "GCF_000001635.27"
    assert row["repository_accession"] == "GEO:GSE168901"
    assert row["geo_id"] == "GSE168901"
    assert row["pubmed_id"] == "25855195"
    assert row["device"] == "OBI:0002630" and row["assay"] == "OBI:0002631"
    assert row["log_transformation"] == "log2 transformation"
    assert row["primary_analysis_completed"] == 1
    assert row["file_format"] == "h5ad"
    assert row["metadata_version"] == "2"
    assert row["metadata_source"] == "uploader"


def test_non_geo_repository_has_no_geo_id(tmp_path, cursor):
    row = save(tmp_path, cursor, repository_accession="ArrayExpress:E-MTAB-11634")
    assert row["repository_accession"] == "ArrayExpress:E-MTAB-11634"
    assert row["geo_id"] is None


def test_v1_fields_still_saved(tmp_path, cursor):
    """Uploads staged with the v1 form (annotation_source / annotation_release_number) still save."""
    row = save(tmp_path, cursor, annotation_source="Ensembl", annotation_release_number="93")
    assert (row["annotation_source"], row["annotation_release"]) == ("Ensembl", "93")


def test_root_term_is_a_curation_flag(tmp_path, cursor):
    row = save(tmp_path, cursor, device="OBI:0000968")
    assert "device_root_term" in json.loads(row["curation_flags"])


def test_ontology_terms_and_keywords_are_saved(tmp_path, cursor):
    cursor.fetchall.return_value = []
    save(tmp_path, cursor, keywords="Atoh1 knockout, hair cells",
         ontology_terms={"tissue_type": ["UBERON:0002227"], "cell_type": ["CL:0000601", "CL:0011113"]})
    inserted = [c.args[1] for c in cursor.execute.call_args_list if "INSERT INTO dataset_ontology_term" in c.args[0]]
    assert inserted == [("DS1", "tissue_type", "UBERON:0002227"),
                        ("DS1", "cell_type", "CL:0000601"), ("DS1", "cell_type", "CL:0011113")]
    tags = [c.args[1] for c in cursor.execute.call_args_list if "INSERT INTO tag " in c.args[0]]
    assert tags == [("Atoh1 knockout",), ("hair cells",)]


def test_pi_name_takes_known_spelling(tmp_path, cursor):
    """A PI name another dataset already has, typed differently, is saved with that dataset's spelling."""
    cursor.fetchone.return_value = ("Jane Doe",)
    row = save(tmp_path, cursor, pi_name=" jane   doe ")
    assert row["pi_name"] == "Jane Doe"
    lookup, = [c.args for c in cursor.execute.call_args_list if c.args[0].startswith("SELECT pi_name")]
    assert lookup[1] == ("jane doe", "DS1")


def test_new_pi_name_is_kept(tmp_path, cursor):
    cursor.fetchone.return_value = None
    assert save(tmp_path, cursor, pi_name="Jane Doe")["pi_name"] == "Jane Doe"


def test_blank_pi_name_is_null(tmp_path, cursor):
    row = save(tmp_path, cursor, pi_name="  ")
    assert row["pi_name"] is None
    assert not any(c.args[0].startswith("SELECT pi_name") for c in cursor.execute.call_args_list)


# ---------------------------------------------------------------- update_dataset_metadata


def test_update_changes_only_given_fields(cursor):
    cursor.fetchone.return_value = (1, None)
    result = md.update_dataset_metadata("DS1", {"title": "New title", "pubmed_id": "PMID:123456"}, source="curator")
    assert result == {"errors": {}, "curation_flags": {}}

    (sql, params), = [c.args for c in cursor.execute.call_args_list if c.args[0].startswith("UPDATE dataset")]
    columns = [a.split(" = ")[0] for a in sql.removeprefix("UPDATE dataset SET ").split(" WHERE")[0].split(", ")]
    assert dict(zip(columns, params)) == {"title": "New title", "pubmed_id": "123456",
                                          "metadata_source": "curator", "curation_flags": None}
    assert params[-1] == "DS1"


def test_update_rejects_invalid_values_without_writing(cursor):
    result = md.update_dataset_metadata("DS1", {"device": "HiSeq", "cell_type": ["hair cell"], "bogus": 1}, source="curator")
    assert set(result["errors"]) == {"device", "cell_type", "bogus"}
    assert not any(c.args[0].startswith("UPDATE") for c in cursor.execute.call_args_list)


def test_update_keeps_flags_for_fields_not_changed(cursor):
    cursor.fetchone.return_value = (1, json.dumps({"annotation_not_loaded": "old", "device_root_term": "old"}))
    result = md.update_dataset_metadata("DS1", {"device": "OBI:0002630"}, source="backfill-llm")
    assert result["curation_flags"] == {"annotation_not_loaded": "old"}


def test_update_replaces_ontology_terms(cursor):
    cursor.fetchone.return_value = (1, None)
    md.update_dataset_metadata("DS1", {"cell_type": ["CL:0000601"]}, source="curator")
    statements = [c.args for c in cursor.execute.call_args_list]
    assert ("DELETE FROM dataset_ontology_term WHERE dataset_id = %s AND field = %s", ("DS1", "cell_type")) in statements
    assert any("INSERT INTO dataset_ontology_term" in sql and params == ("DS1", "cell_type", "CL:0000601")
               for sql, params in statements)


def test_update_pi_name_takes_known_spelling(cursor):
    cursor.fetchone.side_effect = [(1, None), ("Jane Doe",)]
    md.update_dataset_metadata("DS1", {"pi_name": "jane doe"}, source="backfill-llm")
    (sql, params), = [c.args for c in cursor.execute.call_args_list if c.args[0].startswith("UPDATE dataset")]
    columns = [a.split(" = ")[0] for a in sql.removeprefix("UPDATE dataset SET ").split(" WHERE")[0].split(", ")]
    assert dict(zip(columns, params))["pi_name"] == "Jane Doe"

@pytest.mark.parametrize("label, expected", [
    ("single-nucleus-rnaseq", "single-nucleus-rnaseq"),
    ("single-cell-rnaseq", "single-cell-rnaseq"),
])
def test_dataset_type_single_nucleus_is_kept_apart(label, expected):
    """Regression for #1345: single-nucleus labels are no longer folded into single-cell-rnaseq."""
    df = pd.DataFrame({"value": [label]}, index=["dataset_type"])
    assert md.get_value_from_df(df, "dataset_type") == expected


def test_single_nucleus_dtype_is_saved(tmp_path, cursor):
    """Regression for #1345: an upload typed as single-nucleus is saved with that dtype."""
    row = save(tmp_path, cursor, dataset_type="single-nucleus-rnaseq")
    assert row["dtype"] == "single-nucleus-rnaseq"
