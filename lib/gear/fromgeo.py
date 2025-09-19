import json

import pandas as pd
import requests


class FromGeo:
    """
    Things this class handles:
        1. read dataset_metadata.xlsx template
        2. validate information
        3. add GEO metadata to existing user's metadata
    """

    @classmethod
    def get_geo_data(cls, geo_id: str | None=None) -> list:
        """
        Using geo_id (GSExxxxx or GSMxxxx), retrieves GEO metadata for that ID.

        Input
        -----
            geo_id = GSE65633
                - This can also be a sample id (GSM1602228).
        Ouput
        -----
            list of metadata fields and fields. Here is an example for GSE65633:
            ['^SERIES = GSE65633',
            '!Series_title = RNA-seq analysis of neonatal mouse cochlear hair cells',
            '!Series_geo_accession = GSE65633',
            '!Series_status = Public on Mar 01 2015',
            '!Series_submission_date = Feb 04 2015',
            '!Series_last_update_date = Apr 27 2018',
            '!Series_pubmed_id = 25855195',
            '!Series_summary = This study examined transcripts that are enriched in neonatal mouse cochlear hair cells. Hair cells were purified by FACS sorting for GFP fluorescence from the cochleas of transgenic mice in which the endogenous Atoh1 gene was fused with GFP',
            '!Series_overall_design = Two replicates of GFP+ hair cells were compared with all other cochlear cell types that were GFP-',
            '!Series_type = Expression profiling by high throughput sequencing',
            ...
            '!Series_platform_id = GPL13112',
            '!Series_platform_organism = Mus musculus',
            '!Series_platform_taxid = 10090',
            '!Series_sample_organism = Mus musculus',
            '!Series_sample_taxid = 10090',
            '!Series_relation = BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA274564',
            '!Series_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRP053189',
            '']

        Notes
        -----
        Link with GEO data example: https://www.ncbi.nlm.nih.gov/geo/info/download.html
            targ=self - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65633&targ=self&view=full&form=text
            targ=all - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65633&targ=all&view=full&form=text
            targ=gsm - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65633&targ=gsm&view=full&form=text
        """
        if geo_id is None:
            raise Exception("No 'geo_id' provided. Provide GEO GSExxxxx or GSMxxxx ID to continue.")

        #Get metadata from GEO
        r =  requests.get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+str(geo_id)+'&targ=self&view=full&form=text')
        if r.status_code != 200:
            raise Exception("Not able to get GEO metadata. Status code: " + str(r.status_code))

        #Convert GEO data into list
        content = r.content.decode().split('\r\n')
        return content

    @classmethod
    def process_geo_data(cls, content: list | None=None, json_or_dataframe: str | None=None) -> str | pd.DataFrame | None:
        """
        Use this to process the request GET content retrieved from in method get_geo_data()

        This method outputs 2 data structures depending on the setting of argument 'json_or_dataframe'.
            json_or_dataframe='json' - Outputs a JSON formatted string.
            json_or_dataframe='dataframe' - Outputs a pandas dataframe
        Input
        -----
            content = list of content from GEO request
            json_or_dataframe = 'json' or 'dataframe'. Specifies the output format

        Output
        ------
            parsed GEO content is a JSON string or pandas dataframe

        """
        if content is None:
            raise Exception("No 'content' provided. Provide list of metadata from GEO to continue.")
        if json_or_dataframe is None:
            raise Exception("No 'json_or_dataframe' provided. Please specify the output format: 'json' or 'dataframe'.")

        geo_data = {}
        duplicates = {}

        for item in content:
            if item.startswith('!'):
                #Removes !Series, !Sample
                item = item.split('_', 1)[1]

                # Every once in a while there are lines in the file like this:
                #   !sample_table_begin
                if len(item.split(' = ')) > 1:
                    key, value = item.split(' = ', 1)

                    if key not in geo_data:
                        geo_data[key] = value
                    else:
                        if key not in duplicates:
                            duplicates[key] = [value]
                        else:
                            duplicates[key].append(value)

        for k, v in duplicates.items():
            last_value = geo_data[k]
            v.append(last_value)

            # geo_data[k] = v
            geo_data[k] = ", ".join([str(val) for val in v])

        if json_or_dataframe.lower() == 'json':
            return json.dumps(geo_data)

        if json_or_dataframe.lower() == 'dataframe':
            return pd.DataFrame.from_dict(geo_data, orient='index')

    @classmethod
    def add_geo_data(cls, metadata: pd.DataFrame | None=None, geo_data: pd.DataFrame | None=None) -> pd.DataFrame | None:
        """
        If a value is empty (np.nan) in the user's metadata, populate it with the
        value from the same field in the GEO dataframe.

        Input:
            metadata = metadata object (pandas DataFrame)
            geo_data = GEO data object (pandas DataFrame)

        Output: metadata object containing newly added metadata from GEO data object

        NOTE: Does not overwrite fields from original metadata object
        """

        if metadata is None:
            raise Exception("No 'metadata' provided. Provide pandas dataframe containing GEO metadata to continue.")
        if geo_data is None:
            raise Exception("No 'geo_data' provided. Provide pandas dataframe containing GEO metadata to continue.")

        # Get fields in metadata that are 1) Empty and 2) populated from the GEO data
        if 'filled_by_geo' in metadata:
            #print("DEBUG: filled_by_geo in metadata", file=sys.stderr)
            missing_values = metadata[metadata['value'].isna() & (metadata['filled_by_geo'].notna())]

            # For each empty field, fill it with the value from the GEO dataframe (if it's present in the GEO dataframe)
            for index, value in missing_values.itertuples():
                if index in geo_data.index:
                    metadata.loc[index, 'value'] = geo_data.loc[index][0]
        else:
            #print("DEBUG: filled_by_geo not in metadata", file=sys.stderr)
            for index, value in geo_data.itertuples():
                #print("\tDEBUG: iterating a row", file=sys.stderr)

                if index in metadata.index and metadata.loc[index, 'value'] == "":
                    #print("\tDEBUG: assigning metadata value for index {0}".format(index), file=sys.stderr)
                    metadata.loc[index, 'value'] = geo_data.loc[index][0]

        #print("DEBUG: returning metadata object", file=sys.stderr)
        return metadata
