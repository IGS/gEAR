"""gear.metadata_fields and MetadataValidator.validate_values: the v2 dataset metadata template (#390)."""

import pytest
from gear import metadata_fields as mf
from gear.metadatavalidator import MetadataValidator

VALID = {
    "contact_name": "Jane Doe",
    "contact_email": "jane@example.org",
    "contact_institute": "IGS",
    "title": "A dataset",
    "summary": "What was done",
    "dataset_type": "single-cell-rnaseq",
    "reference_annot_id": "Ensembl:93",
    "sample_taxid": "NCBI:txid10090",
    "device": "OBI:0002630",
    "library_cons_method": "Kit A",
    "rna_extraction_method": "Kit B",
    "assay": "OBI:0002631",
    "normalization_method": "TPM",
    "log_transformation": "No Transformation",
    "primary_analysis_completed": "No",
}


@pytest.mark.parametrize("name, value, expected", [
    ("sample_taxid", "10090", "NCBI:txid10090"),
    ("sample_taxid", "txid10090", "NCBI:txid10090"),
    ("sample_taxid", " ncbi:txid9606 ", "NCBI:txid9606"),
    ("sample_taxid", 9606, "NCBI:txid9606"),
    ("sample_taxid", "mouse", "mouse"),
    ("pubmed_id", "25855195", "PMID:25855195"),
    ("pubmed_id", "pmid: 25855195", "PMID:25855195"),
    ("device", "obi_0002630", "OBI:0002630"),
    ("reference_annot_id", "ensembl: 93", "Ensembl:93"),
    ("reference_annot_id", "refseq:GCF_000001635.27", "RefSeq:GCF_000001635.27"),
    ("repository_accession", "GSE168901", "GEO:GSE168901"),
    ("repository_accession", "arrayexpress:E-MTAB-11634", "ArrayExpress:E-MTAB-11634"),
    ("dataset_type", "Single cell RNA-Seq", "single-cell-rnaseq"),
    ("dataset_type", "Spatial", "spatial"),
    ("dataset_type", "Epigenome", "gosling"),
    ("dataset_type", "epigenomics", "gosling"),
    ("dataset_type", "gosling", "gosling"),
    ("log_transformation", "log2 TRANSFORMATION", "log2 transformation"),
    ("primary_analysis_completed", "yes", "Yes"),
    ("title", "  padded  ", "padded"),
    ("title", "", None),
    ("title", float("nan"), None),
    ("title", "None", None),
])
def test_normalize_value(name, value, expected):
    assert mf.normalize_value(name, value) == expected


def test_parsers():
    assert mf.parse_taxon_id("NCBI:txid10090") == 10090
    assert mf.parse_taxon_id("mouse") is None
    assert mf.parse_pubmed_id("PMID:25855195") == "25855195"
    assert mf.parse_reference_annotation("Ensembl:93") == ("Ensembl", "93")
    assert mf.parse_reference_annotation("GenBank:GCA_000001635.9") == ("GenBank", "GCA_000001635.9")
    assert mf.parse_reference_annotation(None) == (None, None)
    assert mf.parse_repository_accession("GEO:GSE1") == ("GEO", "GSE1")
    assert mf.parse_yes_no("Yes") == 1 and mf.parse_yes_no("no") == 0 and mf.parse_yes_no(None) is None


def test_split_keywords_dedupes_case_insensitively():
    assert mf.split_keywords("Atoh1 knockout, expression profile ,atoh1 KNOCKOUT,,") == ["Atoh1 knockout", "expression profile"]
    assert mf.split_keywords(None) == []


def test_upgrade_legacy_values():
    """v1 spreadsheet / form fields are mapped onto their v2 names."""
    upgraded = mf.upgrade_legacy_values({
        "annotation_source": "ensembl", "annotation_release_number": 93.0,
        "geo_accession": " GSE12345 ", "tags": "a, b", "taxon_id": "9606",
    })
    assert upgraded["reference_annot_id"] == "Ensembl:93"
    assert upgraded["repository_accession"] == "GEO:GSE12345"
    assert upgraded["keywords"] == "a, b"
    assert upgraded["sample_taxid"] == "9606"


def test_upgrade_legacy_values_keeps_v2_values():
    upgraded = mf.upgrade_legacy_values({"reference_annot_id": "RefSeq:GCF_1", "annotation_source": "Ensembl",
                                         "annotation_release_number": 93})
    assert upgraded["reference_annot_id"] == "RefSeq:GCF_1"


def test_template_dataset_types_use_labels():
    assert mf.TEMPLATE_DATASET_TYPES == ("single-cell-rnaseq", "bulk-rnaseq", "microarray", "spatial", "epigenome")
    assert all(mf.normalize_value("dataset_type", t) in mf.FIELDS_BY_NAME["dataset_type"].choices
               for t in mf.TEMPLATE_DATASET_TYPES)


def test_required_attributes_follow_the_template():
    required = set(MetadataValidator().required_atts)
    assert {"contact_email", "sample_taxid", "title", "summary", "reference_annot_id", "device", "assay"} <= required
    assert "microarray_hybridization_protocol" not in required
    assert "pubmed_id" not in required


def test_valid_metadata_has_no_errors():
    assert MetadataValidator.validate_values(VALID) == ({}, {})


def test_shortcuts_are_accepted():
    errors, _ = MetadataValidator.validate_values({**VALID, "sample_taxid": "10090", "pubmed_id": "123456",
                                                   "repository_accession": "GSE1", "device": "obi:0002630"})
    assert errors == {}


@pytest.mark.parametrize("field, value", [
    ("contact_email", "not-an-email"),
    ("reference_annot_id", "93"),
    ("reference_annot_id", "UCSC:hg38"),
    ("repository_accession", "SRA:SRP1"),
    ("sample_taxid", "mouse"),
    ("pubmed_id", "PMC123"),
    ("device", "Illumina NovaSeq"),
    ("assay", "OBI:123"),
    ("dataset_type", "proteomics"),
    ("log_transformation", "ln"),
    ("primary_analysis_completed", "maybe"),
    ("title", "x" * 256),
])
def test_invalid_values_are_errors(field, value):
    errors, _ = MetadataValidator.validate_values({**VALID, field: value})
    assert set(errors) == {field}


def test_missing_required_fields():
    errors, _ = MetadataValidator.validate_values({})
    assert set(errors) == {f.name for f in mf.FIELDS if f.required}
    assert all(message == "This field is required." for message in errors.values())


def test_microarray_hybridization_protocol_required_only_for_microarray():
    errors, _ = MetadataValidator.validate_values({**VALID, "dataset_type": "microarray"})
    assert set(errors) == {"microarray_hybridization_protocol"}

    errors, _ = MetadataValidator.validate_values({**VALID, "dataset_type": "microarray",
                                                   "microarray_hybridization_protocol": "Affymetrix GeneChip"})
    assert errors == {}


def test_keyword_limits():
    errors, _ = MetadataValidator.validate_values({**VALID, "keywords": "a, b, c, d, e, f, g"})
    assert set(errors) == {"keywords"}
    errors, _ = MetadataValidator.validate_values({**VALID, "keywords": "x" * 56})
    assert set(errors) == {"keywords"}


def test_obi_root_terms_are_warnings_not_errors():
    errors, warnings = MetadataValidator.validate_values({**VALID, "device": mf.OBI_DEVICE_ROOT, "assay": mf.OBI_ASSAY_ROOT})
    assert errors == {}
    assert set(warnings) == {"device", "assay"}
