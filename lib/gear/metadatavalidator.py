"""
metadatavalidator.py - Validation of dataset upload metadata values.

Provides MetadataValidator, whose static methods check fields from the metadata
spreadsheet uploaded alongside an expression data file.
"""

import pandas as pd
import os, sys
import re
import requests

from gear import metadata_fields as mf


class MetadataValidator:
    """
    Contains methods that validate values within the metadata.xlsx that is uploaded
    along with the expression data file
    """
    def __init__( self ):
        self.required_atts = [f.name for f in mf.FIELDS if f.required]

    @staticmethod
    def validate_values(values: dict) -> tuple[dict, dict]:
        """
        Check every field of the v2 metadata template (see gear.metadata_fields).

        Values are normalized first, so shortcuts like a bare taxon ID are accepted.

        Returns (errors, warnings), each a dict of field name -> message. Errors must be fixed
        before the dataset can be saved; warnings are passed on to the user and curators.
        """
        values = {name: mf.normalize_value(name, value) for name, value in values.items()}
        errors = {}
        warnings = {}

        for f in mf.FIELDS:
            value = values.get(f.name)

            if value is None:
                if f.is_required(values):
                    errors[f.name] = "This field is required."
                continue

            if f.max_length is not None and len(value) > f.max_length:
                errors[f.name] = f"Must be {f.max_length} characters or fewer."
            elif f.choices and value not in f.choices:
                errors[f.name] = "Must be one of: {0}.".format(", ".join(f.choices))
            elif f.pattern and not re.fullmatch(f.pattern, value):
                errors[f.name] = f"Expected {f.pattern_hint}."

        keywords = mf.split_keywords(values.get("keywords"))
        if len(keywords) > mf.MAX_KEYWORDS:
            errors["keywords"] = f"Provide at most {mf.MAX_KEYWORDS} keywords."
        elif any(len(k) > mf.MAX_KEYWORD_LENGTH for k in keywords):
            errors["keywords"] = f"Each keyword must be {mf.MAX_KEYWORD_LENGTH} characters or fewer."

        # Parts of these values are stored in their own, narrower columns
        _, release = mf.parse_reference_annotation(values.get("reference_annot_id"))
        if "reference_annot_id" not in errors and release and len(release) > mf.ANNOTATION_RELEASE_MAX_LENGTH:
            errors["reference_annot_id"] = f"The version after the source must be {mf.ANNOTATION_RELEASE_MAX_LENGTH} characters or fewer."
        repository, accession = mf.parse_repository_accession(values.get("repository_accession"))
        if "repository_accession" not in errors and repository == "GEO" and len(accession) > mf.GEO_ID_MAX_LENGTH:
            errors["repository_accession"] = f"GEO accessions must be {mf.GEO_ID_MAX_LENGTH} characters or fewer."

        if values.get("device") == mf.OBI_DEVICE_ROOT:
            warnings["device"] = "The OBI 'device' root term was given. The curator team will contact you to find a more specific term."
        if values.get("assay") == mf.OBI_ASSAY_ROOT:
            warnings["assay"] = "The OBI 'assay' root term was given. The curator team will contact you to find a more specific term."

        return errors, warnings

    # check that required fields are populated
    @staticmethod
    def validate_required_field(value: str | None=None) -> bool:
        """
        Return True if the value is present and longer than one character.
        """
        is_valid = False
        # NaN never compares equal to anything (including np.nan), so check it with pd.isna
        if value is None or (pd.api.types.is_scalar(value) and pd.isna(value)):
            return is_valid
        else:
            if len(str(value)) > 1:
                is_valid = True

        return is_valid

    @staticmethod
    def validate_tags(value: str | None=None) -> bool:
        """
        Validate the optional semicolon-separated tags field; currently always returns True.
        """
        #Tags are optional so empty is okay
        is_valid = True
        if value is None:
            return is_valid
        else:
            if len(str(value)) > 1:
                if ';' in str(value):
                    is_valid = True
        return is_valid

    @staticmethod
    def validate_email(email: str | None=None) -> bool:
        """
        Return True if the email is present and matches a basic email address pattern.
        """
        #Check the format of the email
        is_valid = False
        if email is None:
            return is_valid
        else:
            if len(str(email)) > 1:
                # Is email correctly formatted? regex from http://emailregex.com/
                if re.match(r"(^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$)", str(email)):
                    is_valid = True

        return is_valid


    # check if pubmed id is valid through URL search
    @staticmethod
    def validate_pubmed_id(pubmed_id: str | None=None) -> bool:
        """
        Return True if the PubMed ID is absent (optional) or numeric.
        """
        is_valid = False
        print("DEBUG: Going to validate this pubmed ID:({0})".format(pubmed_id), file=sys.stderr)
        if pubmed_id is None:
            is_valid = True
        elif str(pubmed_id).isnumeric():
            is_valid = True
        return is_valid


    # check if geo id is valid
    @staticmethod
    def validate_geo_id(geo_id: str | None=None) -> bool:
        """
        Return True if the ID is a GSE accession that resolves on the NCBI GEO website.
        """
        is_valid = False
        if geo_id is None:
            return is_valid

        geo_id = str(geo_id).lower()
        if geo_id.startswith('gse'):
            r=requests.get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=' + geo_id)
            if r.status_code == 200:
                is_valid = True

        return is_valid

    @staticmethod
    def validate_taxon_id(txid: str | None=None) -> bool:
        """
        Currently only checks that the taxon ID is numeric.
        """
        is_valid = False
        if txid is None:
            return is_valid

        if re.match(r"^\d+$", txid):
            is_valid = True
        else:
            hold_txid = re.sub('[^0-9]','', txid)
            if re.match(r"^\d+$", hold_txid):
                is_valid = True

        return is_valid
