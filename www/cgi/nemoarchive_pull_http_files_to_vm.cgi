#!/opt/bin/python3

# nemoarchive_pull_http_files_to_vm.cgi - Transfer an archived dataset from a NeMO Archive HTTPS location into the NeMO Analytics VM
# After the files are pulled, the filenames are extracted and returned to the client

import cgi
import json
import os,sys
from pathlib import Path
from urllib.parse import urlparse

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

import geardb

www_path = Path(__file__).resolve().parents[1]

# Parse gEAR config
# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]
from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
UPLOAD_BASE_DIR = abs_path_www.joinpath("uploads/files")

# TODO: Will need to expand to restricted datasets eventually also

DB_STEP = "pulled_to_vm_status"    # step name in database

def extract_filenames(filename, dest_dir):
    # untar the file and save the filenames to result
    if filename.endswith(".tar.gz"):
        import tarfile
        with tarfile.open(filename, "r:gz") as tar:
            tar.extractall(path=dest_dir)
            return tar.getnames()
    elif filename.endswith(".zip"):
        import zipfile
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall(dest_dir)
            return zip_ref.namelist()
    else:
        return [filename]

def transfer_file(http_path, dest_path):
    # download the file
    import requests
    r = requests.get(http_path, stream=True)
    r.raise_for_status()
    with open(dest_path, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    return dest_path

def pull_http_files_to_vm(http_path, dataset_id):
    s_dataset = geardb.get_submission_dataset_by_dataset_id(dataset_id)
    if not s_dataset:
        raise ValueError("Submission dataset for dataset id {} not found in database".format(dataset_id))

    s_dataset.save_change(attribute=DB_STEP, value="loading")

    parsed_url = urlparse(http_path)
    if parsed_url.scheme != "https":
        raise ValueError("Only HTTPS URLs are supported")
    source_path = parsed_url.path

    dest_dir = Path(source_path).joinpath(dataset_id)
    dest_dir.mkdir(exist_ok=True)
    dest_filename = Path(source_path).name
    result = {"success": False, "filenames":[]}
    try:
        filename = transfer_file(http_path, str(dest_dir.joinpath(dest_filename)))
        result["filenames"] = extract_filenames(filename, dest_dir)
        result["success"] = True
        # Update status in dataset
        s_dataset.save_change(attribute=DB_STEP, value="completed")
    except Exception as e:
        s_dataset.save_change(attribute="log_message", value=str(e))
        s_dataset.save_change(attribute=DB_STEP, value="failed")
        s_dataset.update_downstream_steps_to_canceled(attribute=DB_STEP)
        # NOTE: Original files are not deleted from the "upload" area, so we can try again.
        print(str(e), file=sys.stderr)
        result["success"] = False

    return result

def main():
    form = cgi.FieldStorage()
    http_path = form.getvalue('http_path')
    dataset_id = form.getvalue('dataset_id')

    result = pull_http_files_to_vm(http_path, dataset_id)

    print('Content-Type: application/json\n\n', flush=True)
    print(json.dumps(result))

if __name__ == '__main__':
    main()