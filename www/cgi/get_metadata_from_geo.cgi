#!/opt/bin/python3

"""
Given a passed GEO ID, this returns a JSON object with the metadata for the series.

Returns an empty JSON object if the passed ID wasn't found
"""

import cgi
import os, sys
import json
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

from gear.fromgeo import FromGeo


def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    geo_id = form.getvalue('geo_id')
    
    result = { 'success': 0, 'message': '', 'data': None }

    if geo_id.startswith('GSE'):
        series_content = FromGeo.get_geo_data(geo_id=geo_id)
        series_content = FromGeo.process_geo_data(content=series_content, json_or_dataframe='json')
        series_json = json.loads(series_content)

        # Can be like: GSM1602229, GSM1602230, GSM1602231, GSM1602228
        sample_id_str = series_json['sample_id']

        if ',' in sample_id_str:
            sample_id_list = sample_id_str.split(',')
            first_sample_id = sample_id_list[0]
        elif sample_id_str.startswith('GSM'):
            first_sample_id = sample_id_str

        sample_content = FromGeo.get_geo_data(geo_id=first_sample_id)
        sample_content = FromGeo.process_geo_data(content=sample_content, json_or_dataframe='json')
        sample_json = json.loads(sample_content)

        # Merge the sample metadata in but keep the series accession
        series_json.update(sample_json)
        series_json['geo_accession'] = geo_id
        result['data'] = series_json

    elif geo_id.startswith('GSM'):
        sample_content = FromGeo.get_geo_data(geo_id=geo_id)
        sample_content = FromGeo.process_geo_data(content=sample_content, json_or_dataframe='json')
        sample_json = json.loads(sample_content)
        result['data'] = sample_json

    print(json.dumps(result['data']))

if __name__ == '__main__':
    main()
