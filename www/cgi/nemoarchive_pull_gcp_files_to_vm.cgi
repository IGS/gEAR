#!/opt/bin/python3

# nemoarchive_pull_gcp_files_to_vm.cgi - Run gsutils to pull datasets from a NeMO Archive GCP bucket into the NeMO Analytics VM

import cgi
import json
import os,sys
from pathlib import Path
from google.cloud import storage
import google.auth

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

import geardb

www_path = os.path.abspath(os.path.join('..'))

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
UPLOAD_BASE_DIR = abs_path_www.joinpath("uploads/files")

# TODO: Will need to expand to restricted datasets eventually also
# This bucket should be publicly-accessible without need of specialized credentials
# ! Do not include "gs://"
BUCKET_NAME = "nemo-public"

DB_STEP = "pulled_to_vm_status"    # step name in database

def download_blob(bucket_name, source_blob_name, destination_file_name):
    """Downloads a blob from the bucket."""

    success_dict = {"success":False, "message":"", "filename":""}
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

def pull_gcp_files_to_vm(bucket_path, dataset_id):
    s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
    s_dataset.save_change(attribute=DB_STEP, value="loading")

    source_blob_name = bucket_path.rpartition(BUCKET_NAME + "/")[2] # Retrieve all after the bucket name
    dest_dir = Path(UPLOAD_BASE_DIR).joinpath(dataset_id)
    dest_dir.mkdir(exist_ok=True)
    dest_filename = Path(source_blob_name).name
    result = download_blob(BUCKET_NAME, source_blob_name, str(dest_dir.joinpath(dest_filename)))

    # Update status in dataset
    try:
        if result["success"]:
            s_dataset.save_change(attribute=DB_STEP, value="completed")
    except Exception as e:
        s_dataset.save_change(attribute=DB_STEP, value="failed")
        s_dataset.update_downstream_steps_to_cancelled(attribute=DB_STEP)
        # NOTE: Original files are not deleted from the "upload" area, so we can try again.
        print(str(e), file=sys.stderr)
        result["message"] = "Submission {} - Dataset - {} -- Could not save status to database.".format("test", dataset_id)
        result["success"] = False
    return result

def main():
    form = cgi.FieldStorage()
    bucket_path = form.getvalue('bucket_path')
    dataset_id = form.getvalue('dataset_id')

    result = pull_gcp_files_to_vm(bucket_path, dataset_id)

    print('Content-Type: application/json\n\n', flush=True)
    print(json.dumps(result))

if __name__ == '__main__':
    main()