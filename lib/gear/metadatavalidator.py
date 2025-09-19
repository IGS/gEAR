import pandas as pd
import numpy as np
import os, sys
import re
import requests


class MetadataValidator:
    """
    Contains methods that validate values within the metadata.xlsx that is uploaded
    along with the expression data file
    """
    def __init__( self ):
        self.required_atts = [
            'annotation_release_number',
            'annotation_source',
            'contact_email',
            'contact_institute',
            'contact_name',
            'dataset_type',
            'sample_taxid',
            'summary',
            'title'
        ]

    # check that required fields are populated
    @staticmethod
    def validate_required_field(value: str | None=None) -> bool:
        is_valid = False
        if value is None:
            return is_valid
        else:
            if len(str(value)) > 1 and value != np.nan:
                is_valid = True

        return is_valid

    @staticmethod
    def validate_tags(value: str | None=None) -> bool:
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
