#!/opt/bin/python3

# pull_nemoarchive_gcp_files_to_vm.cgi - Run gsutils to pull datasets from a NeMO Archive GCP bucket into the NeMO Analytics VM

import cgi
import json
import os,sys
from pathlib import Path
from google.cloud import storage
import google.auth

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

www_path = os.path.abspath(os.path.join('..'))

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
UPLOAD_BASE_DIR = abs_path_www.joinpath("uploads/files")
UPLOAD_BASE_DIR = "/tmp"

# TODO: Will need to expand to restricted datasets eventually also
# This bucket should be publicly-accessible without need of specialized credentials
# ! Do not include "gs://"
BUCKET_NAME = "nemo-public"

def download_blob(bucket_name, source_blob_name, destination_file_name):
    """Downloads a blob from the bucket."""

    success_dict = {"success":"", "message":"", "filename":""}
    try:
        # https://google-auth.readthedocs.io/en/latest/user-guide.html
        #TODO maybe put in .env file instead
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = str(Path(www_path).joinpath(this.servercfg["nemoanalytics_import"]["credentials_json"]))
        credentials, project = google.auth.default()
        storage_client = storage.Client(project=project, credentials=credentials)

        bucket = storage_client.bucket(bucket_name, user_project=project)

        # Construct a client side representation of a blob.
        # Note `Bucket.blob` differs from `Bucket.get_blob` as it doesn't retrieve
        # any content from Google Cloud Storage. As we don't need additional data,
        # using `Bucket.blob` is preferred here.
        blob = bucket.blob(source_blob_name)
        blob.download_to_filename(destination_file_name)
        success_dict["success"] = True
        success_dict["message"] = "Downloaded storage object {} from bucket {} to local file {}.".format(
            source_blob_name, bucket_name, destination_file_name
        )
        success_dict["filename"] = destination_file_name
    except Exception as e:
        success_dict["success"] = False
        success_dict["message"] = str(e)

    return success_dict


def main():
    form = cgi.FieldStorage()
    bucket_path = form.getvalue('bucket_path')
    dataset_id = form.getvalue('dataset_id')

    source_blob_name = bucket_path.rpartition(BUCKET_NAME + "/")[2] # Retrieve all after the bucket name
    dest_dir = Path(UPLOAD_BASE_DIR).joinpath(dataset_id)
    dest_dir.mkdir(exist_ok=True)
    dest_filename = Path(source_blob_name).name
    result = download_blob(BUCKET_NAME, source_blob_name, str(dest_dir.joinpath(dest_filename)))
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()