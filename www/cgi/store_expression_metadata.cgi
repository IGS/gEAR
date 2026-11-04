#!/opt/bin/python3

"""
Used by the expression uploader, this stores metadata from the form and saves
it to a file until the expression data are ready to be saved.

Writes a file at: ../uploads/files/<session_id>/<share_uid>/metadata.json
"""

import cgi
import json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from gear import metadata_fields as mf
from gear.metadata import find_curation_flags
from gear.metadatavalidator import MetadataValidator
from werkzeug.utils import secure_filename

def main():
    print('Content-Type: application/json\n\n')
    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    share_uid = form.getfirst('share_uid')

    # Both are used as path components, so refuse (rather than silently rewrite) anything but a plain name
    safe_session_id = secure_filename(session_id or '')
    if not safe_session_id or safe_session_id != session_id:
        print(json.dumps({'success': 0, 'error': 'Invalid session_id'}))
        return

    safe_share_uid = secure_filename(share_uid or '')
    if not safe_share_uid or safe_share_uid != share_uid:
        print(json.dumps({'success': 0, 'error': 'Invalid share_uid'}))
        return

    user = geardb.get_user_from_session_id(session_id)
    if not user:
        print(json.dumps({'success': 0, 'error': 'Invalid session_id'}))
        return

    # Paths are built from the secure_filename() results (identical to the inputs after the checks above)
    user_upload_file_base = "../uploads/files/{0}".format(safe_session_id)

    # make sure the directory exists
    if not os.path.exists(user_upload_file_base):
        os.makedirs(user_upload_file_base)

    user_upload_file_base = os.path.join(user_upload_file_base, safe_share_uid)

    # make sure the directory exists
    if not os.path.exists(user_upload_file_base):
        os.makedirs(user_upload_file_base)

    result = {'success': 0, 'error': None, 'errors': {}, 'warnings': {}}

    # Fields of the metadata template (see lib/gear/metadata_fields.py), stored in their canonical form
    values = {f.name: mf.normalize_value(f.name, form.getfirst(f.name)) for f in mf.FIELDS}

    errors, _warnings = MetadataValidator.validate_values(values)
    if errors:
        result['errors'] = errors
        print(json.dumps(result))
        return

    taxon_id = mf.parse_taxon_id(values['sample_taxid'])
    if geardb.get_organism_id_by_taxon_id(taxon_id) is None:
        result['errors'] = {'sample_taxid': "gEAR doesn't support this organism yet. Please contact the gEAR team to request it."}
        print(json.dumps(result))
        return

    # Things the gEAR team follows up on (#421). The upload continues, and curators are e-mailed
    #  when it is finalized.
    result['warnings'] = find_curation_flags(values)

    formdata = {
            **values,
            'owner_id': user.id,
            'dataset_uid': form.getfirst('dataset_uid'),
            'share_uid': form.getfirst('share_uid'),
            'user_pii_affirmed': form.getfirst('user_pii_affirmed'),
            'metadata_version': mf.METADATA_VERSION,
            'default_plot_type': None,
            'schematic_image': None
    }

    # Add some extra parameters that we can access in the future in the UI
    formdata["perform_primary_analysis"] = False
    formdata["dataset_format"] = ""

    # Save the metadata to a file
    metadata_filename = os.path.join(user_upload_file_base, 'metadata.json')
    try:
        with open(metadata_filename, 'w') as f:
            f.write(json.dumps(formdata))
        result['success'] = 1
    except Exception as e:
        result['error'] = str(e)
        pass

    print(json.dumps(result))

if __name__ == '__main__':
    main()
