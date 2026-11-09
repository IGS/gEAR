"""
metadata.py - Dataset metadata reading, validation, and persistence.

Reads dataset metadata (e.g. from the metadata spreadsheet template), optionally
enriches it with GEO data, validates it, and saves it to MySQL or JSON.
"""

import json
import re
import sys
from io import StringIO

import geardb
import mysql.connector
import numpy as np
import pandas as pd
from gear import metadata_fields as mf
from gear.fromgeo import FromGeo
from gear.metadatavalidator import MetadataValidator as mdv

# Fields of the dataset_ontology_term table that are rolled up from sample-level metadata (#422)
ONTOLOGY_TERM_FIELDS = ("cell_type", "tissue_type")

# Curation flags (see find_curation_flags) and the metadata field each one is about
CURATION_FLAG_FIELDS = {
    "annotation_not_loaded": "reference_annot_id",
    "device_root_term": "device",
    "assay_root_term": "assay",
}


def is_na(value):
    """
    An empty cell in a pandas dataframe returns as NaN, which is a numpy float type.
    Simple check df.loc(index_field, column) is None does not evaluate.
    If I change an empty cell (NaN) to a string. the value now becomes 'nan', which I can
    evaluate against.
    """
    if str(value) == 'nan' or value == np.nan or len(str(value)) < 1:
        return True
    else:
        return False

def get_value_from_df(df, index_label):
    """
    Given a dataframe and an index row label, return the value from column 'value',
    or return None if no value is found.

    Note: Want to return None for the MySQL insert
    """
    if index_label not in df.index:
        #raise Exception("KeyError: label '" + index_label +"' not found in index." )
        return None

    value = df.loc[index_label, 'value']

    if isinstance(value, dict):
        if 'value' in value:
            value = value['value']

    # Clean up dtype labels
    if index_label == 'dataset_type':
        if 'single' in str(value).lower():
            value = 'single-cell-rnaseq'
        if 'bulk' in str(value).lower():
            value = 'bulk-rnaseq'
        if value == 'microarray':
            value = 'microarray'
        if value == 'ChIP-Seq':
            value = 'chip-seq'
        if value == 'ATAC-Seq':
            value = 'atac-seq'

    # If value is NaN, return None
    if is_na(value) is True:
        return None
    else:
        return value



class Metadata:
    """
    Things this class handles:
        1. read dataset metadata file (XLS or JSON)
        2. validate information
        3. upload metadata to gEAR MySQL
    """
    def __init__(self, metadata=None, file_path: str | None=None):
        self.metadata = metadata

        if file_path is None:
            raise Exception("No 'file_path' provided. Please provide the path to your metadata file to continue.")

        self.file_path = file_path
        self.read_file()

    def make_spatial_h5ad_adjustment(self, dataset_format):
        """
        Adjusts the 'dtype' field in the metadata to 'spatial-h5ad' if the dataset format is not 'spatial'
        but the current metadata dtype is 'spatial'. This is to accommodate spatial datasets that are
        uploaded in h5ad format, so that they are not routed to the SpatialAnalysis class

        Parameters
        ----------
        dataset_format : str
            The format of the dataset to check against the metadata dtype.

        Raises
        ------
        Exception
            If self.metadata is None, indicating that no metadata is available.
        """
        if self.metadata is None:
            raise Exception("No metadata found in self.metadata. Provide JSON or Excel metadata template file to continue.")
        dtype = self.metadata.loc['dataset_type', 'value']

        if dataset_format != "spatial" and dtype == "spatial":
            self.metadata.loc['dataset_type', 'value'] = 'spatial-h5ad'

    def read_file(self) -> None:
        """
        Reads dataset_metadata.xlsx or dataset_metadata.json into a pandas dataframe

        Input
        -----
            filepath - /path/to/your/metadata.xlsx or
                       /path/to/your/metadata.json

        Output
        ------
            Populates the metadata attribute of the object as a pandas
            dataframe indexed on 'fields' and has column 'values'
        """
        file_path = self.file_path

        #Read in file
        if file_path.endswith('xlsx') or file_path.endswith('xls'):
            try:
                self.metadata = pd.read_excel(file_path, sheet_name='metadata', index_col='field')
            except Exception as err:
                raise Exception(err)
        elif file_path.endswith('json'):
            try:
                json_data = {'field': [], 'value': []}

                with open(file_path) as json_file:
                    data = json.loads(json_file.read())

                    for d in data:
                        json_data['field'].append(d)
                        json_data['value'].append(data[d])

                pd_df = pd.read_json(StringIO(json.dumps(json_data)), orient='columns')
                self.metadata = pd_df.set_index('field')

            except Exception as err:
                raise Exception(err)
        else:
            raise Exception("Unrecognized metadata file extension in file: {}".format(file_path))

        self._upgrade_legacy_fields()

    def _upgrade_legacy_fields(self) -> None:
        """
        Add v2 template fields filled from their v1 equivalents (see
        metadata_fields.upgrade_legacy_values), so v1 spreadsheets and uploads staged before
        the v2 template was deployed can still be read. Also turns a template dataset_type label
        into the stored value.
        """
        if self.metadata is None:
            return

        current = self.values()
        upgraded = mf.upgrade_legacy_values(current)
        for field, value in upgraded.items():
            if current.get(field) is None and value is not None:
                self.metadata.loc[field, 'value'] = value

        # The template shows some types by label (e.g. epigenome); the uploader form needs the stored value
        dataset_type = mf.normalize_value('dataset_type', current.get('dataset_type'))
        if dataset_type is not None:
            self.metadata.loc['dataset_type', 'value'] = dataset_type

    def values(self) -> dict:
        """Return the metadata as a dict of field -> value, with blank values as None."""
        if self.metadata is None:
            return {}
        return {field: get_value_from_df(self.metadata, field) for field in self.metadata.index}


    def populate_from_geo(self):
        """
        populates metadata fields where 'filled_by_geo' values are 'Y'. If a GEO GSExxxx ID is not
        given, then GEO metadata will not be retrieved.

        Note: This only populates fields that are empty (loaded as NaN by pandas).
              User values are not overwritten.

        Input
        -----
            self = MetadataUploader object where self.metadata is a pandas dataframe populated from
                    an user's metadata template file.

        Ouput
        -----
            Sets self.metadata = pandas DataFrame where empty fields are now populated by GEO info
        """
        if self.metadata is None:
            raise Exception("No metadata found in self.metadata. Provide JSON or Excel metadata template file to continue.")

        prefix, geo_series_id = mf.parse_repository_accession(get_value_from_df(self.metadata, 'repository_accession'))

        if prefix != 'GEO' or not geo_series_id:
            #Stop here if a GEO ID is not given
            return self

        # Get series metadata from GEO
        series_content = FromGeo.get_geo_data(geo_id=str(geo_series_id))

        # Convert data into pandas dataframe
        try:
            series_df = FromGeo.process_geo_data(content=series_content, json_or_dataframe='dataframe')
            if not isinstance(series_df, pd.DataFrame):
                raise Exception(f"Expected DataFrame from process_geo_data, got {type(series_df)}")
        except Exception as err:
            raise Exception("Unable to process GEO data for series {0}. Error: {1}".format(geo_series_id, err))

        # The v2 template has no 'filled_by_geo' column; mark the fields GEO may fill instead
        if 'filled_by_geo' not in self.metadata:
            self.metadata['filled_by_geo'] = [
                'Y' if field in mf.FIELDS_BY_NAME and mf.FIELDS_BY_NAME[field].filled_by_geo else np.nan
                for field in self.metadata.index
            ]

        # Combine user's metdata and series metadata
        updated_metadata = FromGeo.add_geo_data(metadata=self.metadata, geo_data=series_df)

        # Repeat to get sample metadata
        # sample_id = series_df.loc['sample_id', 0][0]
        sample_ids = series_df.loc['sample_id', 0]
        sample_id = str(sample_ids).split(',', 1)[0]

        samp_content = FromGeo.get_geo_data(geo_id=sample_id)
        samp_df = FromGeo.process_geo_data(content=samp_content, json_or_dataframe='dataframe')
        if not isinstance(samp_df, pd.DataFrame):
            raise Exception(f"Expected DataFrame from process_geo_data, got {type(samp_df)}")
        self.metadata = FromGeo.add_geo_data(metadata=self.metadata, geo_data=samp_df)


    def validate(self):
        """
        Runs validation checks on the metadata values (see MetadataValidator.validate_values).

        This actually updates the metadata DataFrame with error messages for
        any invalid fields, and marks required fields in the 'is_required' column.

        returns True if passed all validity tests
        """

        if self.metadata is None:
            raise Exception("No values to evaluate. Please load a metadata file first.")

        values = self.values()
        errors, _warnings = mdv.validate_values(values)

        for f in mf.FIELDS:
            if f.is_required(values) and f.name in self.metadata.index:
                self.metadata.loc[f.name, 'is_required'] = 1

        for field, message in errors.items():
            self.metadata.loc[field, 'message'] = message
            print(f"Invalid metadata field {field}: {message}", file=sys.stderr)

        return not errors


    def add_field_value(self, field=None, value=None):
        """
        Intended to adding the gEAR MySQL fields like dataset_uid, share_uid, owner_id
        to the metadata dataframe object
        """
        if self.metadata is None:
            raise Exception("No values to evaluate. Please load a metadata file first.")
        if field is None:
            raise Exception("No 'field' given. Please supply a field.")
        # if value is None:
        #     raise Exception("No 'value' given. Please supply a value.")

        # Add field and value to metadata
        self.metadata.loc[field, 'value'] = value

        return self

    def save_to_mysql(self, status=None, is_public=0, source='uploader') -> None:
        """
        Saves metadata to gEAR MySQL table 'dataset'. If present, also saves keywords (as tags)
        and the ontology terms rolled up from the sample-level metadata.

        Notes:
            is_public = 0
                All datasets will save as private. Once the upload is complete,
                the user can change the dataset to public on the dataset manager.
            load_status = 'completed'
                All datasets will save as 'completed'.
            source
                Stored in dataset.metadata_source to record where the metadata came from.

        Afterwards self.curation_flags holds any issues a curator should follow up on
        (see find_curation_flags).
        """
        if self.metadata is None:
            raise Exception("No values to evaluate. Please load a metadata file first.")

        if status is None:
            status = 'completed'

        df = self.metadata
        values = self.values()

        cnx = geardb.Connection()
        cursor = cnx.get_cursor()

        dataset_uid = get_value_from_df(df, 'dataset_uid')
        owner_id = str( get_value_from_df(df, 'owner_id') )
        dataset_title = get_value_from_df(df, 'title')

        # Get organism gEAR ID using taxon_id
        organism_taxid = mf.parse_taxon_id(get_value_from_df(df, 'sample_taxid'))
        if organism_taxid is None:
            raise Exception("No taxon_id found in metadata. Please provide a taxon_id.")
        organism_id = geardb.get_organism_id_by_taxon_id(organism_taxid)
        if organism_id is None:
            raise Exception("No organism found with taxon_id: {0}. Please add the organism to gEAR before uploading this dataset.".format(organism_taxid))

        # Store a missing GEO ID as NULL, not as "" or "None" (#1106)
        repository_prefix, repository_id = mf.parse_repository_accession(get_value_from_df(df, 'repository_accession'))
        repository_accession = mf.normalize_value('repository_accession', get_value_from_df(df, 'repository_accession'))
        geo_id = repository_id if repository_prefix == 'GEO' else None

        ldesc = get_value_from_df(df, 'summary')
        dtype = get_value_from_df(df, 'dataset_type')
        schematic_image = get_value_from_df(df, 'schematic_image')
        share_uid = get_value_from_df(df, 'share_uid')
        default_data_format = 'raw'
        has_h5ad = 1

        if dtype == "spatial":
            has_h5ad = 0

        pubmed_id = get_value_from_df(df, 'pubmed_id')
        if isinstance(pubmed_id, str) and pubmed_id.strip() in ('', 'None'):
            pubmed_id = None

        if isinstance(pubmed_id, str):
            pubmed_id = re.sub(r'PMID:?\s*', '', pubmed_id.strip(), flags=re.IGNORECASE)

            # Users entering multiple pubmed IDs will cause failure.  Take the first
            # one instead and append the rest to the Long description.
            pubmed_id = pubmed_id.replace(' ', ',')
            pubmed_ids = pubmed_id.split(',')
            pubmed_ids = [i for i in pubmed_ids if len(i) > 4]

            if len(pubmed_ids):
                pubmed_id = pubmed_ids.pop(0)
            else:
                pubmed_id = None

            if len(pubmed_ids):
                if not isinstance(ldesc, str):
                    ldesc = str(ldesc) if ldesc is not None else ''

                ldesc += "<br>Additional Pubmed IDS: {0}".format(', '.join(pubmed_ids))

        # Legacy (v1 template) fields, still filled in by GEO lookups or older spreadsheets
        platform_id = get_value_from_df(df, 'platform_id')
        instrument_model = get_value_from_df(df, 'instrument_model')
        library_selection = get_value_from_df(df, 'library_selection')
        library_source = get_value_from_df(df, 'library_source')
        library_strategy = get_value_from_df(df, 'library_strategy')

        contact_email = get_value_from_df(df, 'contact_email')
        contact_institute = get_value_from_df(df, 'contact_institute')
        contact_name = get_value_from_df(df, 'contact_name')
        pi_name = match_known_pi_name(cursor, mf.normalize_value('pi_name', get_value_from_df(df, 'pi_name')), dataset_uid)
        annotation_source, annotation_release = mf.parse_reference_annotation(get_value_from_df(df, 'reference_annot_id'))
        default_plot_type = get_value_from_df(df, 'default_plot_type')
        user_pii_affirmed = get_value_from_df(df, 'user_pii_affirmed')

        if user_pii_affirmed is None:
            user_pii_affirmed = 0

        # v2 template fields
        device = mf.normalize_value('device', get_value_from_df(df, 'device'))
        assay = mf.normalize_value('assay', get_value_from_df(df, 'assay'))
        library_cons_method = get_value_from_df(df, 'library_cons_method')
        rna_extraction_method = get_value_from_df(df, 'rna_extraction_method')
        microarray_hybridization_protocol = get_value_from_df(df, 'microarray_hybridization_protocol')
        normalization_method = get_value_from_df(df, 'normalization_method')
        log_transformation = mf.normalize_value('log_transformation', get_value_from_df(df, 'log_transformation'))
        primary_analysis_completed = mf.parse_yes_no(get_value_from_df(df, 'primary_analysis_completed'))
        file_format = get_value_from_df(df, 'dataset_format')

        self.curation_flags = find_curation_flags(values, organism_id=organism_id)
        curation_flags = json.dumps(self.curation_flags) if self.curation_flags else None

        add_dataset_sql = """
        REPLACE INTO dataset (id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, ldesc, date_added, dtype, schematic_image, share_id, math_default, load_status, has_h5ad, platform_id, instrument_model, library_selection, library_source, library_strategy, contact_email, contact_institute, contact_name, annotation_source, annotation_release, plot_default, user_pii_affirmed, repository_accession, device, assay, library_cons_method, rna_extraction_method, microarray_hybridization_protocol, normalization_method, log_transformation, primary_analysis_completed, file_format, metadata_version, metadata_source, curation_flags, pi_name)
        VALUES              (%s, %s,       %s,    %s,          %s,        %s,     %s,        %s,    NOW(),      %s,    %s,              %s,       %s,           %s,          %s,       %s,          %s,               %s,                %s,             %s,               %s,            %s,                %s,           %s,                %s,                 %s,           %s,                %s,                   %s,     %s,    %s,                  %s,                    %s,                                %s,                   %s,                 %s,                         %s,          %s,               %s,              %s,             %s)
        """

        # Insert dataset info to database
        try:
            cursor.execute(add_dataset_sql, (dataset_uid, owner_id, dataset_title, organism_id, pubmed_id,
                                             geo_id, is_public, ldesc, dtype, schematic_image,
                                             share_uid, default_data_format, status, has_h5ad, platform_id,
                                             instrument_model, library_selection, library_source, library_strategy, contact_email,
                                             contact_institute, contact_name, annotation_source, annotation_release, default_plot_type,
                                             user_pii_affirmed, repository_accession, device, assay, library_cons_method,
                                             rna_extraction_method, microarray_hybridization_protocol, normalization_method, log_transformation, primary_analysis_completed,
                                             file_format, mf.METADATA_VERSION, source, curation_flags, pi_name,))
            cnx.commit()
        except mysql.connector.Error as err:
            raise Exception("Failed to insert metadata: {0}".format(err))

        # Keywords were called 'tags' in the v1 template, and are stored as tags
        save_dataset_keywords(cursor, cnx, dataset_uid, mf.split_keywords(get_value_from_df(df, 'keywords')))

        ontology_terms = get_value_from_df(df, 'ontology_terms')
        if isinstance(ontology_terms, dict):
            save_dataset_ontology_terms(cursor, cnx, dataset_uid, ontology_terms)

        cursor.close()


    def write_json(self, file_path=None):
        """
        Writes the metadata DataFrame to a JSON file or returns it as a JSON string.

        If no file path is provided, non-metadata columns ('value', 'message', 'is_required') are removed before exporting.
        If a file path is provided, only the 'value' column is retained for export.

        Args:
            file_path (str, optional): The path to the output JSON file. If None, returns the JSON string.

        Returns:
            str: The JSON representation of the metadata if file_path is None; otherwise, writes to the specified file.

        Raises:
            Exception: If metadata is not loaded before calling this method.
        """

        if self.metadata is None:
            raise Exception("No values to evaluate. Please load a metadata file first.")

        for col in self.metadata.columns.to_numpy():
            if file_path is None:
                # Remove non-metadata columns
                if col != 'value' and col != 'message' and col != 'is_required':
                    self.metadata = self.metadata.drop(col, axis=1)
            else:
                # Only keep metadata fields and values
                if col != 'value':
                    self.metadata = self.metadata.drop(col, axis=1)

        # Write metadata to json file
        return self.metadata.to_json(path_or_buf=file_path, orient='index')


def is_reference_annotation_loaded(source: str | None, release: str | None, organism_id: int | None) -> bool:
    """
    Return True if gEAR has genes loaded for the given reference annotation and organism.

    Only Ensembl releases are tracked (gene.ensembl_release). gEAR doesn't record which
    GenBank/RefSeq assembly genes came from, so those can't be confirmed and return False.
    """
    if source != 'Ensembl' or organism_id is None or not str(release or '').isdigit():
        return False

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    try:
        cursor.execute("SELECT 1 FROM gene WHERE organism_id = %s AND ensembl_release = %s LIMIT 1",
                       (organism_id, int(release)))
        return cursor.fetchone() is not None
    finally:
        cursor.close()
        cnx.close()


def find_curation_flags(values: dict, organism_id: int | None = None) -> dict:
    """
    Look for submitted metadata a gEAR curator should follow up on (#421).

    Returns a dict of flag -> message for the user. Flags are:
        annotation_not_loaded   reference annotation isn't (confirmed to be) loaded in gEAR
        device_root_term        the OBI 'device' root term was given instead of a specific device
        assay_root_term         the OBI 'assay' root term was given instead of a specific assay
    """
    flags = {}

    if organism_id is None:
        taxon_id = mf.parse_taxon_id(values.get('sample_taxid'))
        if taxon_id is not None:
            organism_id = geardb.get_organism_id_by_taxon_id(taxon_id)

    source, release = mf.parse_reference_annotation(values.get('reference_annot_id'))
    if source and not is_reference_annotation_loaded(source, release, organism_id):
        flags['annotation_not_loaded'] = ("The reference annotation {0}:{1} is not loaded in gEAR yet. You can continue "
                                          "with the upload, and the gEAR team will look into adding it.").format(source, release)

    _errors, warnings = mdv.validate_values(values)
    for field in ('device', 'assay'):
        if field in warnings:
            flags[f'{field}_root_term'] = warnings[field]

    return flags


def match_known_pi_name(cursor, pi_name: str | None, dataset_id: str) -> str | None:
    """
    Return the spelling another dataset already uses for this PI name, so "jane doe" is saved as
    "Jane Doe" if that's how it was entered before. Otherwise return the name as given.

    The match ignores case because the dataset table's collation does. dataset_id is left out of the
    search, so a dataset's own old spelling doesn't win over a correction to it.
    """
    if not pi_name:
        return pi_name
    cursor.execute("SELECT pi_name FROM dataset WHERE pi_name = %s AND id != %s LIMIT 1", (pi_name, dataset_id))
    row = cursor.fetchone()
    return row[0] if row and row[0] else pi_name


def save_dataset_keywords(cursor, cnx, dataset_id: str, keywords: list[str]) -> None:
    """Link keywords to a dataset as tags, adding any tags that don't exist yet."""
    if not keywords:
        return

    #Get list of tags already in gEAR
    cursor.execute("SELECT label, id FROM tag")
    cached_tags = {row[0].lower(): row[1] for row in (cursor.fetchall() or [])}

    add_tag_sql = """
        INSERT INTO tag (label)
        VALUES (%s)
    """
    add_datasettag_sql = """
        INSERT INTO dataset_tag (tag_id, dataset_id)
        VALUES (%s, %s)
    """
    for keyword in keywords:
        #Check if new tag is already in database
        if keyword.lower() not in cached_tags:
            #New tag. Add it to database and keep its id
            cursor.execute(add_tag_sql, (keyword,))
            cnx.commit()
            tag_id = cursor.lastrowid
        else:
            #Tag exists. Get its id
            tag_id = cached_tags[keyword.lower()]

        cursor.execute(add_datasettag_sql, (tag_id, dataset_id,) )
        cnx.commit()


def save_dataset_ontology_terms(cursor, cnx, dataset_id: str, terms_by_field: dict) -> None:
    """
    Replace a dataset's rows in dataset_ontology_term for each field given.

    terms_by_field is like {'cell_type': ['CL:0000601'], 'tissue_type': ['UBERON:0002227']}.
    Fields not given are left alone.
    """
    for field, term_ids in terms_by_field.items():
        if field not in ONTOLOGY_TERM_FIELDS:
            raise ValueError(f"Unknown ontology term field: {field}")

        cursor.execute("DELETE FROM dataset_ontology_term WHERE dataset_id = %s AND field = %s", (dataset_id, field))
        for term_id in sorted(set(term_ids or [])):
            cursor.execute("INSERT INTO dataset_ontology_term (dataset_id, field, term_id) VALUES (%s, %s, %s)",
                           (dataset_id, field, term_id))
    cnx.commit()


def update_dataset_metadata(dataset_id: str, values: dict, source: str) -> dict:
    """
    Update metadata fields of an existing dataset.

    This is the entry point for editing metadata after upload, including backfilling older
    datasets with the v2 template fields (e.g. after a curator accepts suggested values).
    Only the fields given are changed. Field names are those of the v2 template
    (gear.metadata_fields), plus 'cell_type' / 'tissue_type' as lists of ontology term IDs.

    source is recorded in dataset.metadata_source, e.g. 'curator' or 'backfill-llm'.

    Returns {'errors': {field: message}, 'curation_flags': {flag: message}}. Nothing is saved
    if there are errors.
    """
    values = mf.upgrade_legacy_values(values)
    errors, _warnings = mdv.validate_values(values)
    # Fields that weren't given aren't being changed, so they can't be missing
    errors = {field: message for field, message in errors.items() if field in values}

    for field in ONTOLOGY_TERM_FIELDS:
        for term_id in values.get(field) or []:
            if not re.fullmatch(r"[A-Za-z]+:\d+", str(term_id)):
                errors[field] = f"Not an ontology term ID: {term_id}"

    unknown = set(values) - set(mf.FIELDS_BY_NAME) - set(ONTOLOGY_TERM_FIELDS) - mf.LEGACY_FIELDS
    for field in unknown:
        errors[field] = "Unknown metadata field."

    if errors:
        return {'errors': errors, 'curation_flags': {}}

    columns = {}
    for name, value in values.items():
        if name not in mf.FIELDS_BY_NAME:
            continue
        value = mf.normalize_value(name, value)

        if name == 'sample_taxid':
            organism_id = geardb.get_organism_id_by_taxon_id(mf.parse_taxon_id(value))
            if organism_id is None:
                return {'errors': {name: f"No organism in gEAR has taxon ID {value}."}, 'curation_flags': {}}
            columns['organism_id'] = organism_id
        elif name == 'reference_annot_id':
            columns['annotation_source'], columns['annotation_release'] = mf.parse_reference_annotation(value)
        elif name == 'repository_accession':
            prefix, accession = mf.parse_repository_accession(value)
            columns['repository_accession'] = value
            columns['geo_id'] = accession if prefix == 'GEO' else None
        elif name == 'pubmed_id':
            columns['pubmed_id'] = mf.parse_pubmed_id(value)
        elif name == 'primary_analysis_completed':
            columns['primary_analysis_completed'] = mf.parse_yes_no(value)
        elif mf.FIELDS_BY_NAME[name].column:
            columns[mf.FIELDS_BY_NAME[name].column] = value

    columns['metadata_source'] = source

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    try:
        cursor.execute("SELECT organism_id, curation_flags FROM dataset WHERE id = %s", (dataset_id,))
        row = cursor.fetchone()
        if row is None:
            return {'errors': {'dataset_id': f"No dataset with ID {dataset_id}."}, 'curation_flags': {}}
        organism_id, existing_flags = row
        curation_flags = find_curation_flags(values, organism_id=columns.get('organism_id', organism_id))

        # Re-check only the flags whose fields are being changed, and keep the others
        existing_flags = json.loads(existing_flags) if existing_flags else {}
        for flag, field in CURATION_FLAG_FIELDS.items():
            if field in values:
                existing_flags.pop(flag, None)
        curation_flags = {**existing_flags, **curation_flags}
        columns['curation_flags'] = json.dumps(curation_flags) if curation_flags else None

        if 'pi_name' in columns:
            columns['pi_name'] = match_known_pi_name(cursor, columns['pi_name'], dataset_id)

        # Column names come from gear.metadata_fields, never from the caller
        assignments = ", ".join(f"{column} = %s" for column in columns)
        cursor.execute(f"UPDATE dataset SET {assignments} WHERE id = %s", (*columns.values(), dataset_id))

        if 'keywords' in values:
            cursor.execute("DELETE FROM dataset_tag WHERE dataset_id = %s", (dataset_id,))
            save_dataset_keywords(cursor, cnx, dataset_id, mf.split_keywords(values['keywords']))

        terms = {field: values[field] for field in ONTOLOGY_TERM_FIELDS if field in values}
        if terms:
            save_dataset_ontology_terms(cursor, cnx, dataset_id, terms)

        cnx.commit()
    finally:
        cursor.close()
        cnx.close()

    return {'errors': {}, 'curation_flags': curation_flags}
