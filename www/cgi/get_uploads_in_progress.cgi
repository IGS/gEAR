#!/opt/bin/python3

"""
For a given session_id, returns information on any uploads 
in an incomplete state.

Data structure returned:

{
   'success': 1,
   'uploads': [
        { 
   ]
}

"""

import cgi, json
import os, sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb



def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')

    result = {'success':0, 'uploads':[], 'message':''}

    if not session_id:
        result['message'] = 'No session_id provided'
        print(json.dumps(result))
        return
    
    user_upload_file_base = "../uploads/files/{0}".format(session_id)

    # If this directory doesn't exist, there are no uploads in progress
    if not os.path.exists(user_upload_file_base):
        result['success'] = 1
        print(json.dumps(result))
        return
    
    # Get the list of share_dirs in the directory
    #  Each of these should be a share ID for an upload in progress
    share_ids = os.listdir(user_upload_file_base)
    for share_id in share_ids:
        share_dir = "{0}/{1}".format(user_upload_file_base, share_id)
        metadata_file = "{0}/metadata.json".format(share_dir)
        
        # get some attributes from the metadata file
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

            result['uploads'].append( {
                    'share_id': share_id,
                    'dataset_id': metadata.get('dataset_uid', ''),
                    'dataset_type': metadata.get('dataset_type', ''),
                    'title': metadata.get('title', ''),
                    'status': 'metadata uploaded',
                    'load_step': 'upload-dataset'
                }
            )

        tarball_data_file = "{0}/{1}.tar.gz".format(share_dir, share_id)

        if os.path.exists(tarball_data_file):
            result['uploads'][-1]['status'] = 'datafile uploaded'
            result['uploads'][-1]['load_step'] = 'process-dataset'

        processing_status_json_file = os.path.join(share_dir, 'status.json')

        if os.path.exists(processing_status_json_file):
            with open(processing_status_json_file, 'r') as f:
                status_json = json.load(f)
                processing_status = status_json.get('status', '')

                if processing_status == 'processing':
                    result['uploads'][-1]['status'] = 'processing'
                    result['uploads'][-1]['load_step'] = 'process-dataset'
                elif processing_status == 'complete':
                    result['uploads'][-1]['status'] = 'processed'
                    result['uploads'][-1]['load_step'] = 'finalize-dataset'
    
    result['success'] = 1
    print(json.dumps(result))


if __name__ == '__main__':
    main()
