import ast
import json
import numpy as np
import pandas as pd
import os, sys, re
import requests
import mysql.connector

import geardb
from gear.fromgeo import FromGeo
from gear.metadatavalidator import MetadataValidator as mdv


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
        if value == 'single-cell RNA-Seq': value = 'single-cell-rnaseq'
        if value == 'bulk RNA-Seq': value = 'bulk-rnaseq'
        if value == 'microarray': value = 'microarray'
        if value == 'ChIP-Seq': value = 'chip-seq'
        if value == 'ATAC-Seq': value = 'atac-seq'

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
    def __init__(self, metadata=None, file_path=None):
        self.metadata = metadata
        self.file_path = file_path

        if self.file_path is not None:
            self.read_file(file_path=file_path)


    def read_file(self, file_path=None):
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
                    data = ast.literal_eval(json_file.read())
                    for d in data:
                        json_data['field'].append(d)
                        json_data['value'].append(data[d])

                pd_df = pd.read_json(json.dumps(json_data), orient='columns')
                self.metadata = pd_df.set_index('field')

            except Exception as err:
                raise Exception(err)
        else:
            raise Exception("Unrecognized metadata file extension in file: {}".format(file_path))


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
            raise Exception("No metadata found in self.metadata. Provide read an excel metadata template file to continue.")

        geo_series_id = self.metadata.loc['geo_accession', 'value']

        if isinstance(geo_series_id, dict):
            if 'value' in geo_series_id:
                geo_series_id = geo_series_id['value'].strip()
            else:
                raise Exception("geo_series_id was a dict but didn't have 'value' as a key")
        else:
            if isinstance(geo_series_id, str):
                geo_series_id = geo_series_id.strip()

        if is_na(geo_series_id) is True:
            #Stop here if geo_id is not given
            return self

        # Get series metadata from GEO
        series_content = FromGeo.get_geo_data(geo_id=geo_series_id)

        # Convert data into pandas dataframe
        series_df = FromGeo.process_geo_data(content=series_content, json_or_dataframe='dataframe')

        # Combine user's metdata and series metadata
        updated_metadata = FromGeo.add_geo_data(metadata=self.metadata, geo_data=series_df)

        # Repeat to get sample metadata
        # sample_id = series_df.loc['sample_id', 0][0]
        sample_ids = series_df.loc['sample_id', 0]
        sample_id = sample_ids.split(',', 1)[0]

        samp_content = FromGeo.get_geo_data(geo_id=sample_id)
        samp_df = FromGeo.process_geo_data(content=samp_content, json_or_dataframe='dataframe')
        self.metadata = FromGeo.add_geo_data(metadata=self.metadata, geo_data=samp_df)


    def validate(self):
        """
        Runs validation checks on the metadata values. Checks are imported from
        metadatavalidator.py containing class MetadataValidator and methods

        This actually updates the metadata DataFrame with error messages for
        any invalid fields.

        returns True if passed all validity tests
        """

        if self.metadata is None:
            raise Exception("No values to evaluate. Please load a metadata file first.")

        df = self.metadata

        validator = mdv()
        is_valid = True

        #check for empty required fields
        for idx, v in df.iterrows():
            if idx in validator.required_atts:
                df.loc[idx, 'is_required'] = 1

                if is_na(v['value']):
                    is_valid = False
                    df.loc[idx, 'message'] = "This field is required. "
                    print("This field is required - {}".format(idx), file=sys.stderr)

        #check email format
        #email = df.loc['contact_email', 'value']
        #if is_na(email) is True:
        #    df.loc['contact_email', 'message'] = "This field is required. Please provide a valid email. "
        #else:
        #    is_email_valid = mdv.validate_email(email)
        #    if is_email_valid is False:
        #        is_valid = False
        #        df.loc['contact_email', 'message'] = "This field is required. The email does not seem valid. Please check it is typed correctly. "

        #check pubmed_id
        #pubmed_id = df.loc['pubmed_id', 'value']
        #if is_na(pubmed_id) is False:
        #    is_pubmed_valid = mdv.validate_pubmed_id(pubmed_id)
        #    if is_pubmed_valid is False:
        #       is_valid = False
        #        df.loc['pubmed_id', 'message'] = "Unable to confirm with URL search. Please ensure the value entered is correct. "

        #check geo_id
        #geo_id = df.loc['geo_accession', 'value']
        #if is_na(geo_id) is False:
        #    is_geo_valid = mdv.validate_geo_id(geo_id)
        #    if is_geo_valid is False:
        #        is_valid = False
        #        df.loc['geo_accession', 'message'] = "Unable to confirm with URL search. Please ensure the value entered is correct. "

        # check the taxon ID
        #taxon_id = str(df.loc['taxon_id', 'value'])
        if 'sample_taxid' in df.index:
            taxon_id = str(df.loc['sample_taxid', 'value'])
        elif 'taxon_id' in df.index:
            taxon_id = str(df.loc['taxon_id', 'value'])
        else:
            raise Exception("No taxon id found")
        is_taxon_valid = mdv.validate_taxon_id(taxon_id)
        if is_taxon_valid is False:
            is_valid = False
            df.loc['sample_taxid', 'message'] = "Taxon ID should be numeric (only). Please ensure the value entered is correct. "

        self.metadata = df
        return is_valid


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

    def get_field_value(self, field=None):
        """
        Accessor for field attributes in the metadata dataframe.
        """
        fv = self.metadata.loc[field, 'value']
        if isinstance(fv, dict):
            if 'value' in fv:
                fv = fv['value']
            else:
                raise Exception("Field value is a dict() but has no 'value' key")

        return fv

    def save_to_mysql(self, status=None):
        """
        Saves metadata to gEAR MySQL table 'dataset'. If present, also saves tags.

        Notes:
            is_public = 0
                All datasets will save as private. Once the upload is complete,
                the user can change the dataset to public on the dataset manager.
            load_status = 'pending'
                All datasets will save as 'pending'.
        """
        if self.metadata is None:
            raise Exception("No values to evaluate. Please load a metadata file first.")

        if status is None:
            status = 'pending'

        df = self.metadata

        cnx = geardb.Connection()
        cursor = cnx.get_cursor()

        dataset_uid = get_value_from_df(df, 'dataset_uid')
        owner_id = str( get_value_from_df(df, 'owner_id') )
        dataset_title = get_value_from_df(df, 'title')

        # Get organism gEAR ID using taxon_id
        organism_id = None
        organism_taxid = get_value_from_df(df, 'sample_taxid')
        organism_qry = ( "SELECT id FROM organism WHERE taxon_id = %s" )
        cursor.execute(organism_qry, ( str(organism_taxid), ))
        for (id, ) in cursor:
            organism_id = id

        geo_id = str( get_value_from_df(df, 'geo_accession') ).strip()

        # Make all datasets private if not specified in the metadata
        try:
            is_public = get_value_from_df(df, 'is_public')
        except:
            is_public = 0

        ldesc = get_value_from_df(df, 'summary')
        dtype = get_value_from_df(df, 'dataset_type')
        schematic_image = get_value_from_df(df, 'schematic_image')
        share_uid = get_value_from_df(df, 'share_uid')
        default_data_format = 'raw'
        has_h5ad = 1

        pubmed_id = str( get_value_from_df(df, 'pubmed_id') ).strip()

        # Users entering multiple pubmed IDs will cause failure.  Take the first
        # one instead and append the rest to the Long description.
        pubmed_id = pubmed_id.replace(' ', ',')
        pubmed_ids = pubmed_id.split(',')
        pubmed_ids = [i for i in pubmed_ids if len(i) > 3]
        pubmed_id = pubmed_ids.pop()

        if len(pubmed_ids):
            ldesc += "<br>Additional Pubmed IDS: {0}".format(', '.join(pubmed_ids))

        platform_id = get_value_from_df(df, 'platform_id')
        instrument_model = get_value_from_df(df, 'instrument_model')
        library_selection = get_value_from_df(df, 'library_selection')
        library_source = get_value_from_df(df, 'library_source')
        library_strategy = get_value_from_df(df, 'library_strategy')
        contact_email = get_value_from_df(df, 'contact_email')
        contact_institute = get_value_from_df(df, 'contact_institute')
        contact_name = get_value_from_df(df, 'contact_name')
        annotation_source = get_value_from_df(df, 'annotation_source')
        annotation_release = get_value_from_df(df, 'annotation_release_number')
        default_plot_type = get_value_from_df(df, 'default_plot_type')

        add_dataset_sql = """
        INSERT INTO dataset (id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, ldesc, date_added, dtype, schematic_image, share_id, math_default, load_status, has_h5ad, platform_id, instrument_model, library_selection, library_source, library_strategy, contact_email, contact_institute, contact_name, annotation_source, annotation_release, plot_default)
        VALUES              (%s, %s,       %s,    %s,          %s,        %s,     %s,        %s,    NOW(),      %s,    %s,              %s,       %s,           %s,          %s,       %s,          %s,               %s,                %s,             %s,               %s,            %s,                %s,           %s,                %s,                 %s)
        """

        # Insert dataset info to database
        try:
            cursor.execute(add_dataset_sql, (dataset_uid, owner_id, dataset_title, organism_id, pubmed_id,
                                             geo_id, is_public, ldesc, dtype, schematic_image,
                                             share_uid, default_data_format, status, has_h5ad, platform_id,
                                             instrument_model, library_selection, library_source, library_strategy, contact_email,
                                             contact_institute, contact_name, annotation_source, annotation_release, default_plot_type,))
            cnx.commit()
        except mysql.connector.Error as err:
            raise Exception("Failed to insert metadata: {0}".format(err))

        #Handle any tags user might have included
        tags = get_value_from_df(df, 'tags')

        # if tags is not None or tags != np.nan:
        if is_na(tags) is False and tags is not None:
            tag_list = []

            if isinstance(tags, str):
                if ',' in tags:
                    raw_tags = tags.split(', ')
                else:
                    raw_tags = [tags]
            elif isinstance(tags, list):
                raw_tags = tags

            # Ensures duplicates are removed
            for tag in raw_tags:
                if tag not in tag_list:
                    tag_list.append(tag)

            #Get list of tags already in gEAR
            qry_get_tags = """
                SELECT label, id
                FROM tag
            """
            cached_tags = {}
            cursor.execute(qry_get_tags)
            for row in cursor:
                cached_tags[row[0].lower()] = row[1]

            add_tag_sql = """
                INSERT INTO tag (label)
                VALUES (%s)
            """
            add_datasettag_sql = """
                INSERT INTO dataset_tag (tag_id, dataset_id)
                VALUES (%s, %s)
            """
            for tag in tag_list:
                #Check if new tag is already in database
                if tag.lower() not in cached_tags:
                    #New tag. Add it to database and keep its id
                    cursor.execute(add_tag_sql, (tag,))
                    cnx.commit()
                    tag_id = cursor.lastrowid
                else:
                    #Tag exists. Get its id
                    tag_id = cached_tags[tag.lower()]

                cursor.execute(add_datasettag_sql, (tag_id, dataset_uid,) )
                cnx.commit()

        cursor.close()

    #TODO: This is not currently used anywhere.
    def _list_invalid_fields(self):
        """
        Returns a multi-line string listing the invalid metadata fields

        Note: Perform after validation is completed via method validate().
        """

        if 'is_valid' not in self.metadata.columns:
            raise Exception("No column 'is_valid' in dataframe. Please first run method validate().")

        msg = ''
        false_fields = self.metadata.loc[self.metadata['is_valid'] == False]
        ffcount, cols = false_fields.shape

        if ffcount <= 0:
            msg += "Looks good. No invalid metadata fields found."
        else:
            msg += "The following fields were found to contain invalid values:\n"
            for i, row in false_fields.iterrows():
                msg += "\tField: " + str(i) + " has invalid value: " + str(row[0]) + ".\n"

        return msg


    def write_json(self, file_path=None):
        """
        Writes the metadata fields and corresponding values to JSON.

        If filepath is given, the JSON will write to disk where filepath specifies.
        If filepath is None (not given), the JSON will return string
        """

        if self.metadata is None:
            raise Exception("No values to evaluate. Please load a metadata file first.")


        if file_path is None:
            #Remove non-metadata columns
            for col in self.metadata.columns.values:
                if col != 'value' and col != 'message' and col != 'is_required':
                    self.metadata = self.metadata.drop(col, axis=1)

            # Return metadata as JSON
            return self.metadata.to_json(path_or_buf=file_path, orient='index')
        else:
            #Only keep metadata fields and values
            for col in self.metadata.columns.values:
                if col != 'value':
                    self.metadata = self.metadata.drop(col, axis=1)

            # Write metadata to json file
            self.metadata.to_json(path_or_buf=file_path, orient='index')
            return self
