#!/opt/bin/python3

"""
Used by the expression uploader, this stores metadata from the form and saves
it to a file until the expression data are ready to be saved.

Writes a file at: ../uploads/files/<session_id>_<share_uid>.json
"""

import cgi
import json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')

    user_upload_file_base = '../uploads/files'

    user = geardb.get_user_from_session_id(session_id)
    result = {'success': 0}

    # names are changed here so the files are compatible with the legacy ones
    formdata = {
            'owner_id': user.id,
            'dataset_uid': form.getvalue('dataset_uid'),
            'share_uid': form.getvalue('share_uid'),
            'title': form.getvalue('title'),
            'summary': form.getvalue('summary'),
            'dataset_type': form.getvalue('dataset_type'),
            'annotation_source': form.getvalue('annotation_source'),
            'annotation_release_number': form.getvalue('annotation_version'),
            'geo_accession': form.getvalue('geo_id'),
            'contact_name': form.getvalue('contact_name'),
            'contact_email': form.getvalue('contact_email'),
            'contact_institute': form.getvalue('contact_institute'),
            'sample_taxid': form.getvalue('taxon_id'),
            'sample_organism': form.getvalue('organism'),
            'platform_id': form.getvalue('platform_id'),
            'instrument_model': form.getvalue('instrument'),
            'library_selection': form.getvalue('library_selection'),
            'library_source': form.getvalue('library_source'),
            'library_strategy': form.getvalue('library_strategy'),
            'pubmed_id': form.getvalue('pubmed_id'),
            # These needed to be added/supported for real
            'expression_unit': 'normalized log count',
            'tags': None,
            'default_plot_type': None,
            'schematic_image': None
    }

    # Save the metadata to a file
    metadata_filename = os.path.join(user_upload_file_base, session_id + '_' + form.getvalue('share_uid') + '.json')
    try:
        with open(metadata_filename, 'w') as f:
            f.write(json.dumps(formdata))
        result['success'] = 1
    except Exception as e:
        pass

    print(json.dumps(result))

if __name__ == '__main__':
    main()