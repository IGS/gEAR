"""
metadata_fields.py - Definition of the dataset-level metadata template (v2, GitHub issue #390).

This is the single source of truth for the fields a user fills in when uploading a dataset.
It is used by:
  - gear.metadatavalidator, to check submitted values
  - gear.metadata, to normalize values and map them to columns of the 'dataset' table
  - bin/build_metadata_templates.py, to build the downloadable Excel template and the
    field list the uploader form reads (www/user_templates/metadata_fields.json)

Fields that gEAR fills in itself (cell_type, tissue_type, file_format) are not listed here.
"""

import re
import unicodedata
from dataclasses import asdict, dataclass, field

METADATA_VERSION = "2"

DATASET_TYPES = ("single-cell-rnaseq", "bulk-rnaseq", "microarray")

# Types gEAR supports beyond the FAIR template's list
GEAR_ONLY_DATASET_TYPES = ("spatial", "gosling")

# What the template and uploader show for a stored type, where the two differ. Users know
# Gosling data as epigenome tracks, not by the viewer's name.
DATASET_TYPE_LABELS = {"gosling": "epigenome"}

# Choices offered in the template's dataset_type dropdown
TEMPLATE_DATASET_TYPES = tuple(DATASET_TYPE_LABELS.get(t, t) for t in DATASET_TYPES + GEAR_ONLY_DATASET_TYPES)

LOG_TRANSFORMATIONS = ("No Transformation", "log2 transformation", "log10 transformation")

YES_NO = ("Yes", "No")

REFERENCE_ANNOTATION_SOURCES = ("Ensembl", "GenBank", "RefSeq")

REPOSITORY_PREFIXES = ("GEO", "ArrayExpress", "CNGB", "NeMO")

# Root terms of the OBI 'device' and 'assay' branches. Users who can't find a specific term are
#  told to enter these, and a curator follows up with them.
OBI_DEVICE_ROOT = "OBI:0000968"
OBI_ASSAY_ROOT = "OBI:0000070"

MAX_KEYWORDS = 6

# tag.label is VARCHAR(55)
MAX_KEYWORD_LENGTH = 55

# Some fields are split up, and only part of the value is stored in a dataset column. These are
#  the widths of those columns (see create_schema.sql), so values that won't fit are rejected at the
#  metadata step instead of failing when the dataset is saved (GitHub #1323, #340).
ANNOTATION_RELEASE_MAX_LENGTH = 50  # dataset.annotation_release, the part after "Ensembl:", "RefSeq:", ...
GEO_ID_MAX_LENGTH = 50              # dataset.geo_id, the part of a GEO accession after "GEO:"
PUBMED_ID_MAX_LENGTH = 20           # dataset.pubmed_id, the part after "PMID:"


@dataclass(frozen=True)
class MetadataField:
    """One row of the dataset metadata template."""
    name: str
    instructions: str
    example: str = ""
    required: bool = False
    # Required only when another field has the given value, e.g. {"dataset_type": "microarray"}
    required_if: dict = field(default_factory=dict)
    max_length: int | None = None
    # Regex the (normalized) value must match in full
    pattern: str | None = None
    # Human readable description of the pattern, shown in error messages
    pattern_hint: str = ""
    choices: tuple = ()
    # Column of the 'dataset' table this field is stored in, if stored directly
    column: str | None = None
    # Filled in from GEO (by matching field name) when left blank and a GEO accession is given
    filled_by_geo: bool = False

    def is_required(self, values: dict) -> bool:
        """Return True if this field is required, given the other submitted values."""
        if self.required:
            return True
        return any(values.get(other) == wanted for other, wanted in self.required_if.items())

    def to_json(self) -> dict:
        """Serializable form of this field, for the uploader form."""
        return asdict(self)


FIELDS = (
    MetadataField(
        name="contact_name", filled_by_geo=True,
        instructions="First and last name of the person to contact about the dataset. This can be the principal investigator, or the lab member who prepared the data.",
        required=True, max_length=100, column="contact_name",
    ),
    MetadataField(
        name="contact_email", filled_by_geo=True,
        instructions="Email address of the contact person.",
        required=True, max_length=100, column="contact_email",
        pattern=r"[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+", pattern_hint="a valid email address",
    ),
    MetadataField(
        name="contact_institute", filled_by_geo=True,
        instructions="Name of institution under which the contact person is employed.",
        required=True, max_length=255, column="contact_institute",
    ),
    MetadataField(
        name="pi_name",
        instructions="First and last name of the principal investigator (lab head) responsible for the dataset, used to find a lab's datasets. The uploader form suggests PI names already on this portal; choose one if it matches, so each PI is spelled one way. If the study has more than one PI, give the lead PI. Leave blank if there is no single PI, e.g. for a consortium dataset.",
        example="Jane Doe",
        max_length=100, column="pi_name",
    ),
    MetadataField(
        name="title",
        instructions="A concise title for the dataset. This typically includes the dataset type, the organism studied, and conditions examined. Limit to 255 characters.",
        example="RNA-seq analysis of neonatal mouse cochlear hair cells",
        required=True, max_length=255, column="title",
    ),
    MetadataField(
        name="summary",
        instructions="An abstract describing the goal of the project and the experiments that were performed.",
        example="This study examined transcripts that are enriched in neonatal mouse cochlear hair cells. Hair cells were purified by FACS sorting for GFP fluorescence from the cochleas of transgenic mice in which the endogenous Atoh1 gene was fused with GFP.",
        required=True, max_length=65535, column="ldesc",
    ),
    MetadataField(
        name="dataset_type",
        instructions="Directs data to the correct part of the portal. If your study has more than one type of data, submit each type as a separate dataset. If you are doing single-nucleus RNA-seq or SMART-seq, use single-cell-rnaseq. Use spatial for spatial transcriptomics, and epigenome for epigenome tracks (e.g. ChIP-seq or ATAC-seq).",
        example="single-cell-rnaseq",
        required=True, choices=DATASET_TYPES + GEAR_ONLY_DATASET_TYPES, column="dtype",
    ),
    MetadataField(
        name="reference_annot_id",
        instructions="Reference annotation the data were aligned to, as source:version. Supported sources are Ensembl, GenBank and RefSeq. For Ensembl give the release number; for GenBank and RefSeq give the assembly accession. Examples: Ensembl:93, RefSeq:GCF_000001635.27, GenBank:GCA_000001635.9. If your reference is not yet loaded into this portal you can still upload, and the curator team will contact you about loading it.",
        example="Ensembl:93",
        required=True, max_length=len("GenBank:") + ANNOTATION_RELEASE_MAX_LENGTH,
        pattern=r"(Ensembl|GenBank|RefSeq):\S+", pattern_hint="Ensembl:<release>, GenBank:<accession> or RefSeq:<accession>",
    ),
    MetadataField(
        name="repository_accession",
        instructions="Accession of the data in a public repository, with the repository as prefix: GEO (e.g. GEO:GSE168901), ArrayExpress (e.g. ArrayExpress:E-MTAB-11634), CNGB (e.g. CNGB:CNP0002521) or NeMO. Leave blank if the data have not been submitted to one of these.",
        example="GEO:GSE168901",
        max_length=100, column="repository_accession",
        pattern=r"(GEO|ArrayExpress|CNGB|NeMO):\S+", pattern_hint="GEO:, ArrayExpress:, CNGB: or NeMO: followed by the accession",
    ),
    MetadataField(
        name="sample_taxid",
        instructions="NCBI taxonomy ID of the one species the samples come from, e.g. human NCBI:txid9606, mouse NCBI:txid10090, rat NCBI:txid10116, zebrafish NCBI:txid7955. A dataset can only hold one species. If your study has samples from more than one species, split the data by species and submit each as its own dataset.",
        example="NCBI:txid10090",
        required=True,
        pattern=r"NCBI:txid\d+", pattern_hint="NCBI:txid followed by the numeric taxon ID",
    ),
    MetadataField(
        name="pubmed_id", filled_by_geo=True,
        instructions='PubMed ID of the publication describing the data, with "PMID:" prefix.',
        example="PMID:25855195",
        max_length=len("PMID:") + PUBMED_ID_MAX_LENGTH,
        pattern=r"PMID:\d+", pattern_hint="PMID: followed by the numeric PubMed ID",
    ),
    MetadataField(
        name="keywords",
        instructions=f"Up to {MAX_KEYWORDS} keywords or phrases related to the dataset that are not otherwise captured in the metadata, separated by commas.",
        example="Atoh1 knockout, expression profile",
    ),
    MetadataField(
        name="device",
        instructions=f"Term ID from the Ontology for Biomedical Investigations (OBI) 'device' branch for the device used to generate the primary data, e.g. Illumina NovaSeq 6000 OBI:0002630, Illumina HiSeq 2500 OBI:0002002. If you can't find the term you need, use the branch root term {OBI_DEVICE_ROOT} and the curator team will follow up with you.",
        example="OBI:0002048",
        required=True, column="device",
        pattern=r"OBI:\d{7}", pattern_hint="OBI: followed by the 7-digit term number",
    ),
    MetadataField(
        name="library_cons_method",
        instructions="Name of the kit/protocol used for cDNA library construction or, if a kit is not used, a link to protocols.io or a DOI for the methods description.",
        example="Nextera XT tagmentation and 96-barcode indexing (Illumina)",
        required=True, max_length=65535, column="library_cons_method",
    ),
    MetadataField(
        name="rna_extraction_method",
        instructions="Name of the kit/protocol used for RNA extraction or, if a kit is not used, a link to protocols.io or a DOI for the methods description.",
        example="High Pure RNA Isolation kit (Roche Diagnostics)",
        required=True, max_length=65535, column="rna_extraction_method",
    ),
    MetadataField(
        name="microarray_hybridization_protocol",
        instructions="Name of the kit/protocol used for microarray hybridization or, if a kit is not used, a link to protocols.io or a DOI for the methods description. Required for microarray datasets.",
        example="Affymetrix GeneChip Array",
        required_if={"dataset_type": "microarray"}, max_length=65535, column="microarray_hybridization_protocol",
    ),
    MetadataField(
        name="assay",
        instructions=f"Term ID from the OBI 'assay' branch for the technique used to generate the primary data, e.g. bulk RNA-seq assay OBI:0003090, single-cell RNA sequencing assay OBI:0002631, single-nucleus RNA sequencing assay OBI:0003109, microarray OBI:0400147. If you can't find the term you need, use the branch root term {OBI_ASSAY_ROOT} and the curator team will follow up with you.",
        example="OBI:0002631",
        required=True, column="assay",
        pattern=r"OBI:\d{7}", pattern_hint="OBI: followed by the 7-digit term number",
    ),
    MetadataField(
        name="normalization_method",
        instructions="Normalization method used on the data, e.g. TPM, RPM/CPM, RPKM/FPKM, SCnorm, GeTMM, DESeq2, ComBat-Seq, QN.",
        example="TPM",
        required=True, max_length=100, column="normalization_method",
    ),
    MetadataField(
        name="log_transformation",
        instructions="Transformation applied to the count data.",
        example="No Transformation",
        required=True, choices=LOG_TRANSFORMATIONS, column="log_transformation",
    ),
    MetadataField(
        name="primary_analysis_completed",
        instructions="A primary analysis is the clustering and dimensionality reduction (UMAP, tSNE or PCA) run on the data before uploading, for example in Seurat or Scanpy. Choose Yes if the observation (sample) metadata includes its cluster assignments or coordinates; the portal then checks that it can find them. Otherwise No.",
        example="No",
        required=True, choices=YES_NO,
    ),
)

FIELDS_BY_NAME = {f.name: f for f in FIELDS}


@dataclass(frozen=True)
class SampleField:
    """One column of the sample-level (observation) metadata template (#422)."""
    name: str
    requirement: str
    instructions: str
    example: str = ""


# What to write in a sample metadata cell that has no value, and why: the INSDC missing value
#  vocabulary (https://www.insdc.org/technical-specifications/missing-value-reporting/), plus "unavailable"
#  from the sample metadata template's sex values. Matched ignoring case.
MISSING_VALUES = ("not applicable", "not collected", "not provided", "restricted access", "missing", "unavailable")

# Other ways people write a missing value. These get a warning suggesting MISSING_VALUES instead.
#  ("none" is left out: it's a real value for treatment.)
MISSING_VALUE_STAND_INS = ("n/a", "na", "n.a.", "not available", "unknown", "null", "-", "?")

SEX_VALUES = ("male", "female", "intersex", "indeterminate")

# Columns users provide in their observation metadata. gEAR adds tissue_ontology_name and
#  cell_ontology_name itself, from the term IDs (see gear.ontology).
SAMPLE_FIELDS = (
    SampleField("sample_id", "Required", "Sample identifier; must be unique for each sample (cell or bulk sample).", "sample1"),
    SampleField("tissue_label", "Optional", "Brief label for the tissue, used for plotting. Maximum 20 characters.", "blood"),
    SampleField("tissue_ontology", "Required",
                "Term ID from the Uberon ontology for the tissue, e.g. skin of body UBERON:0002097, "
                "spiral organ of cochlea UBERON:0002227, brainstem UBERON:0002298.", "UBERON:0000178"),
    SampleField("cell_label", "Optional", "Short label for the cell type, used for plotting. Maximum 20 characters.", "PBMC"),
    SampleField("cell_ontology", "Optional",
                "Term ID(s) from the Cell Ontology for the cell type, separated by commas if more than one, e.g. "
                "cochlear outer hair cell CL:0000601, spiral ganglion neuron CL:0011113.", "CL:2000001"),
    SampleField("treatment", "Required",
                "The treatment or condition the sample was exposed to. Provide a brief but descriptive name.", "Control"),
    SampleField("treatment_label", "Optional", "Shortened treatment label for displays. Maximum 20 characters.", "ctr"),
    SampleField("UMAP_1", "Required if primary analysis is completed",
                "Dimensionality reduction coordinates, one column per component. The portal recognizes UMAP_1/UMAP_2, "
                "tSNE_1/tSNE_2 and PCA_1/PCA_2 (or PC_1/PC_2). Coordinates from other methods (e.g. SWNE) can be "
                "included, but must be chosen by hand in the dataset curator to be plotted.", "1.234"),
    SampleField("UMAP_2", "Required if primary analysis is completed", "Second component; see UMAP_1.", "-0.567"),
    SampleField("cell_clusters", "Required if primary analysis is completed and no coordinates are given",
                "Cluster assignment of each sample from your primary analysis. A Seurat seurat_clusters column, "
                "or a scanpy leiden column, is recognized too.", "cluster_3"),
    SampleField("age_value", "Required",
                "Age: a number, a Carnegie stage, or a developmental stage such as P7 or E15.", "P7"),
    SampleField("age_unit", "Required",
                "Unit of a numeric age, e.g. hours, days, months, years, Carnegie, or postnatal stage.", "days"),
    SampleField("age_order", "Required if age is not numeric", "The order in which age values are displayed (1, 2, 3, ...).", "1"),
    SampleField("sex_assigned_at_birth", "Required",
                "One of: " + ", ".join(SEX_VALUES) + ", or a missing value (see the last column).", "female"),
)

SAMPLE_FIELD_NOTE = ("Any number of optional columns can be added, such as current medications, disease at time of "
                     "sampling or genotype. If you use abbreviations in them, include the full name in a separate column.")

MISSING_VALUE_NOTE = ("If a sample has no value for a column, because it doesn't apply or wasn't recorded, enter one "
                      "of: " + ", ".join(MISSING_VALUES) + ". This tells curators the value is missing on purpose. "
                      "Use it in a whole column if your study has no such data, rather than leaving the column out.")


def is_missing_value(value) -> bool:
    """True for the words in MISSING_VALUES, ignoring case and surrounding spaces."""
    return isinstance(value, str) and value.strip().lower() in MISSING_VALUES


def is_missing_value_stand_in(value) -> bool:
    """True for other ways of writing a missing value, like N/A or unknown (see MISSING_VALUE_STAND_INS)."""
    return isinstance(value, str) and value.strip().lower() in MISSING_VALUE_STAND_INS


def is_blank(value) -> bool:
    """True for None, NaN and empty or whitespace-only strings."""
    if value is None:
        return True
    if isinstance(value, float) and value != value:   # NaN
        return True
    return str(value).strip() in ("", "nan", "None")


def holds_only_missing_values(series) -> bool:
    """
    True for an obs column (pandas Series) whose values are all missing value words like "not applicable",
    or other ways of writing them like N/A, apart from blanks. Such a column only says the data is missing,
    so there's nothing in it to plot.
    """
    if hasattr(series, "cat"):
        values = series.cat.categories
    elif series.dtype == object or str(series.dtype) == "string":
        values = series.dropna().unique()
    else:
        return False   # numbers
    values = [v for v in values if not is_blank(v)]
    return bool(values) and all(is_missing_value(v) or is_missing_value_stand_in(v) for v in values)


def remove_invisible_characters(value: str) -> str:
    """
    Remove the invisible characters that come along with text pasted from Word, Excel or web pages
    (GitHub #59): control characters, zero-width spaces and joiners, byte order marks, soft hyphens and
    other formatting marks. Tabs and line breaks are kept, and non-breaking spaces become plain spaces.
    """
    value = value.replace("\u00a0", " ").replace("\u2028", "\n").replace("\u2029", "\n")
    return "".join(ch for ch in value if ch in "\t\n\r" or unicodedata.category(ch) not in ("Cc", "Cf"))


def normalize_value(name: str, value):
    """
    Return a value in its canonical form, accepting the common shortcuts people type
    (bare numbers for taxon and PubMed IDs, missing or differently-cased prefixes).
    Invisible characters are removed, and blank values become None. Values that can't be
    interpreted are returned unchanged, stripped, so validation can report them.
    """
    if is_blank(value):
        return None

    value = remove_invisible_characters(str(value)).strip()
    if not value:
        return None

    if name == "sample_taxid":
        m = re.fullmatch(r"(?:NCBI:)?(?:txid)?\s*(\d+)", value, re.IGNORECASE)
        if m:
            return f"NCBI:txid{m.group(1)}"
    elif name == "pubmed_id":
        m = re.fullmatch(r"(?:PMID:?)?\s*(\d+)", value, re.IGNORECASE)
        if m:
            return f"PMID:{m.group(1)}"
    elif name == "pi_name":
        # Runs of spaces would make the same name look like a different PI
        return re.sub(r"\s+", " ", value)
    elif name in ("device", "assay"):
        m = re.fullmatch(r"OBI[:_]\s*(\d{7})", value, re.IGNORECASE)
        if m:
            return f"OBI:{m.group(1)}"
    elif name == "reference_annot_id":
        return _normalize_prefix(value, REFERENCE_ANNOTATION_SOURCES)
    elif name == "repository_accession":
        # A bare GEO series ID is common enough to accept
        if re.fullmatch(r"GSE\d+", value, re.IGNORECASE):
            return f"GEO:{value.upper()}"
        return _normalize_prefix(value, REPOSITORY_PREFIXES)
    elif name == "dataset_type":
        lowered = value.lower()
        if "single" in lowered:
            return "single-cell-rnaseq"
        if "bulk" in lowered:
            return "bulk-rnaseq"
        if lowered.startswith("epigen"):
            return "gosling"
        if lowered in DATASET_TYPES + GEAR_ONLY_DATASET_TYPES:
            return lowered
    elif name in ("log_transformation", "primary_analysis_completed"):
        for choice in FIELDS_BY_NAME[name].choices:
            if value.lower() == choice.lower():
                return choice

    return value


def _normalize_prefix(value: str, prefixes: tuple) -> str:
    """Fix the case of a known 'Prefix:' and remove spaces around the colon."""
    if ":" not in value:
        return value
    prefix, rest = value.split(":", 1)
    for known in prefixes:
        if prefix.strip().lower() == known.lower():
            return f"{known}:{rest.strip()}"
    return value


def split_keywords(value) -> list[str]:
    """Split the comma-separated keywords field into a de-duplicated list (case-insensitive)."""
    if is_blank(value):
        return []
    if isinstance(value, (list, tuple)):
        items = value
    else:
        items = str(value).split(",")

    keywords = []
    seen = set()
    for item in items:
        item = str(item).strip()
        if item and item.lower() not in seen:
            seen.add(item.lower())
            keywords.append(item)
    return keywords


def parse_taxon_id(value) -> int | None:
    """'NCBI:txid10090' (or '10090', 'txid10090') -> 10090"""
    value = normalize_value("sample_taxid", value)
    if value is None:
        return None
    m = re.fullmatch(r"NCBI:txid(\d+)", value)
    return int(m.group(1)) if m else None


def parse_pubmed_id(value) -> str | None:
    """'PMID:25855195' -> '25855195', the form stored in dataset.pubmed_id"""
    value = normalize_value("pubmed_id", value)
    if value is None:
        return None
    m = re.fullmatch(r"PMID:(\d+)", value)
    return m.group(1) if m else None


def parse_reference_annotation(value) -> tuple[str | None, str | None]:
    """'Ensembl:93' -> ('Ensembl', '93'), stored in dataset.annotation_source / annotation_release"""
    value = normalize_value("reference_annot_id", value)
    if value is None or ":" not in value:
        return None, None
    source, release = value.split(":", 1)
    return source, release


def parse_repository_accession(value) -> tuple[str | None, str | None]:
    """'GEO:GSE168901' -> ('GEO', 'GSE168901')"""
    value = normalize_value("repository_accession", value)
    if value is None or ":" not in value:
        return None, None
    prefix, accession = value.split(":", 1)
    return prefix, accession


def parse_yes_no(value) -> int | None:
    """'Yes' -> 1, 'No' -> 0, anything else -> None"""
    value = normalize_value("primary_analysis_completed", value)
    if value == "Yes":
        return 1
    if value == "No":
        return 0
    return None


# v1 template / v1 uploader form fields that upgrade_legacy_values() reads
LEGACY_FIELDS = {"annotation_source", "annotation_release_number", "annotation_release", "geo_accession", "tags", "taxon_id"}


def upgrade_legacy_values(values: dict) -> dict:
    """
    Return a copy of a metadata dict with fields from the v1 template (and the v1 uploader form)
    mapped onto their v2 names. v2 values already present win. Lets old spreadsheets, and uploads
    that were already in progress when v2 was deployed, still be read.
    """
    values = dict(values)

    def missing(name):
        return is_blank(values.get(name))

    if missing("reference_annot_id") and not is_blank(values.get("annotation_source")):
        release = values.get("annotation_release_number", values.get("annotation_release"))
        if not is_blank(release):
            # Spreadsheet cells may come back as floats, e.g. 93.0
            if isinstance(release, float) and release.is_integer():
                release = int(release)
            source = str(values["annotation_source"]).strip()
            # The v1 form offered 'ncbi', which could be either NCBI resource
            source = {"ensembl": "Ensembl", "genbank": "GenBank", "refseq": "RefSeq"}.get(source.lower(), source)
            values["reference_annot_id"] = f"{source}:{release}"

    if missing("repository_accession") and not missing("geo_accession"):
        values["repository_accession"] = normalize_value("repository_accession", values["geo_accession"])

    if missing("keywords") and not missing("tags"):
        values["keywords"] = values["tags"]

    if missing("sample_taxid") and not missing("taxon_id"):
        values["sample_taxid"] = values["taxon_id"]

    return values

