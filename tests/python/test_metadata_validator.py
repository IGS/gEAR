"""gear.metadatavalidator: validators for fields in the uploaded metadata spreadsheet."""

import numpy as np
import pandas as pd
import pytest
from gear import metadatavalidator
from gear.metadatavalidator import MetadataValidator


def test_required_attributes():
    assert {"contact_email", "sample_taxid", "title", "summary"} <= set(MetadataValidator().required_atts)


@pytest.mark.parametrize("value, expected", [
    (None, False),
    (np.nan, False),
    (float("nan"), False),
    (pd.NA, False),
    (pd.NaT, False),
    ("", False),
    ("x", False),
    ("ab", True),
    ("A normal title", True),
    (9606, True),
    (["a", "b"], True),
])
def test_validate_required_field(value, expected):
    assert MetadataValidator.validate_required_field(value) is expected


@pytest.mark.parametrize("value", [None, "", "x", "one", "one;two"])
def test_validate_tags_is_optional(value):
    assert MetadataValidator.validate_tags(value) is True


@pytest.mark.parametrize("email, expected", [
    ("someone@example.org", True),
    ("first.last+tag@sub-domain.example.co.uk", True),
    ("no-at-sign.example.org", False),
    ("someone@nodot", False),
    ("a", False),
    (None, False),
])
def test_validate_email(email, expected):
    assert MetadataValidator.validate_email(email) is expected


@pytest.mark.parametrize("pubmed_id, expected", [(None, True), ("12345", True), (12345, True), ("PMID12345", False), ("", False)])
def test_validate_pubmed_id(pubmed_id, expected):
    assert MetadataValidator.validate_pubmed_id(pubmed_id) is expected


@pytest.mark.parametrize("taxid, expected", [
    ("9606", True),
    ("txid9606", True),
    ("10090 (mouse)", True),
    ("mouse", False),
    ("", False),
    ("nan", False),
    (None, False),
])
def test_validate_taxon_id(taxid, expected):
    assert MetadataValidator.validate_taxon_id(taxid) is expected


class FakeResponse:
    def __init__(self, status_code):
        self.status_code = status_code


@pytest.mark.parametrize("status, expected", [(200, True), (404, False)])
def test_validate_geo_id_checks_ncbi(monkeypatch, status, expected):
    urls = []
    monkeypatch.setattr(metadatavalidator.requests, "get", lambda url: urls.append(url) or FakeResponse(status))
    assert MetadataValidator.validate_geo_id("GSE12345") is expected
    assert urls == ["https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse12345"]


@pytest.mark.parametrize("geo_id", [None, "GSM12345", "12345"])
def test_validate_geo_id_rejects_non_series_without_a_request(monkeypatch, geo_id):
    monkeypatch.setattr(metadatavalidator.requests, "get", lambda url: pytest.fail("unexpected request"))
    assert MetadataValidator.validate_geo_id(geo_id) is False
