#!/opt/bin/python3

"""
At this point we should have a directory with at least the following:

- JSON file with metadata
- An H5AD file ready to go
- A status.json file

Then, depending on the upload format, we also should have

- If 3tab:
    - A tarball uploaded by the user
- If Excel:
    - An Excel file uploaded by the user

This script does the following:

- Loads the JSON metadata file and stores it in MySQL
- Migrates the H5AD file to the proper directory
- Migrates the original uploaded file so it can be downloaded by users

Returns a status of these steps as JSON data like this:


result = {
    "success": 1,
    "metadata_loaded": 1,
    "h5ad_migrated": 1,
    "userdata_migrated": 1,
    "message": "All steps completed successfully."
}
"""

import cgi
import json
import shutil
import sys
from pathlib import Path

lib_path = Path(__file__).resolve().parents[2] / 'lib'
sys.path.append(str(lib_path))
import geardb
from gear.metadata import Metadata

user_upload_file_path = Path(__file__).resolve().parents[1] / 'uploads' / 'files'
dataset_final_dir = Path(__file__).resolve().parents[1] / 'datasets'

result = {
    "success": 0,
    "metadata_loaded": 0,
    "h5ad_migrated": 0,
    "userdata_migrated": 0,
    "primary_analysis_migrated": 0,
    "message": ""
}

def main() -> dict:
    print('Content-Type: application/json\n\n', flush=True)

    form = cgi.FieldStorage()
    share_uid = form.getvalue('share_uid')
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('dataset_uid')
    dataset_format = form.getvalue('dataset_format')
    dataset_visibility = form.getvalue('dataset_visibility')

    if dataset_format in ["null", None]:
        result['message'] = 'No dataset format provided.'
        return result

    perform_analysis_migration = form.getvalue("perform_analysis_migration", 0)
    # if string, convert to int
    if isinstance(perform_analysis_migration, str):
        perform_analysis_migration = int(perform_analysis_migration)

    if dataset_visibility == 'private':
        is_public = 0
    elif dataset_visibility == 'public':
        is_public = 1
    else:
        result['message'] = 'Invalid dataset visibility.'
        return result

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        result['message'] = 'User ID not found. Please log in to continue.'
        return result

    dataset_upload_dir = user_upload_file_path / session_id / share_uid
    if not dataset_upload_dir.is_dir():
        result['message'] = 'Upload directory not found.'
        return result

    # Load the metadata
    metadata_file = dataset_upload_dir / 'metadata.json'
    if not metadata_file.is_file():
        result['message'] = 'Metadata file not found.'
        return result

    # Load the metadata into the database
    metadata = Metadata(file_path=str(metadata_file))
    try:
        metadata.make_spatial_h5ad_adjustment(dataset_format)
        metadata.save_to_mysql(status='completed', is_public=is_public)
        result['metadata_loaded'] = 1
    except Exception as e:
        result['message'] = 'Error saving metadata to MySQL: {}'.format(str(e))
        return result

    global dataset_final_dir
    # migrate the H5AD file or Zarr store (spatial)
    if dataset_format == "spatial":
        # spatial files go in a subdirectory
        dataset_final_dir = dataset_final_dir / 'spatial'
        zarr_file = dataset_upload_dir / f'{share_uid}.zarr'
        if not zarr_file.is_dir():
            result['message'] = 'Zarr store not found: {}'.format(zarr_file)
            return result

        try:
            shutil.move(zarr_file, dataset_final_dir / f'{dataset_id}.zarr')
            result['h5ad_migrated'] = 1
        except Exception as e:
            result['message'] = 'Error migrating Zarr store: {}'.format(str(e))
            return result
    else:
        h5ad_file = dataset_upload_dir / f'{share_uid}.h5ad'
        if not h5ad_file.is_file():
            result['message'] = 'H5AD file not found: {}'.format(h5ad_file)
            return result
        h5ad_dest = dataset_final_dir / f'{dataset_id}.h5ad'

        try:
            shutil.move(h5ad_file, h5ad_dest)
            result['h5ad_migrated'] = 1
        except Exception as e:
            result['message'] = 'Error migrating H5AD file: {}'.format(str(e))
            return result


    if dataset_format == 'mex_3tab':
        # migrate the tarball
        tarball_file = dataset_upload_dir / f'{share_uid}.tar.gz'
        tarball_dest = dataset_final_dir / f'{dataset_id}.tar.gz'

        #print(f"DEBUG: Attempting to do: mv {tarball_file} {tarball_dest}", file=sys.stderr)
        try:
            shutil.move(tarball_file, tarball_dest)
            result['userdata_migrated'] = 1
        except Exception as e:
            result['message'] = 'Error migrating tarball file: {}'.format(str(e))
            return result

    elif dataset_format == 'excel':
        # migrate the Excel file
        excel_file = dataset_upload_dir / f'{share_uid}.xlsx'
        excel_dest = dataset_final_dir / f'{dataset_id}.xlsx'

        try:
            shutil.move(excel_file, excel_dest)
            result['userdata_migrated'] = 1
        except Exception as e:
            result['message'] = 'Error migrating Excel file: {}'.format(str(e))
            return result

    elif dataset_format == "spatial":
        # migrate the spatial tarball
        spatial_src = dataset_upload_dir / f'{share_uid}.tar.gz'
        spatial_dest = dataset_final_dir / f'{dataset_id}.tar.gz'

        try:
            shutil.move(spatial_src, spatial_dest)
            result['userdata_migrated'] = 1
        except Exception as e:
            result['message'] = 'Error migrating spatial data tarball file: {}'.format(str(e))
            return result
    elif dataset_format == "h5ad":
        # The h5ad is the userdata too
        result['userdata_migrated'] = 1
    else:
        print(f"DEBUG: dataset_format is {dataset_format}", file=sys.stderr)

    # Migrate the primary analysis JSON
    # If the analysis was not created, it is non-fatal
    if perform_analysis_migration == 1:
        analysis_json = dataset_upload_dir / "analysis_pipeline.json"
        if analysis_json.is_file():
            try:
                shutil.move(analysis_json, dataset_final_dir / f"{dataset_id}.pipeline.json")
            except Exception as e:
                result['message'] = 'Error migrating primary analysis JSON: {}'.format(str(e))
                return result

        result['primary_analysis_migrated'] = 1

        # Move preliminary composition plots if they exist. This should not affect the primary analysis migration
        violin_plot = dataset_upload_dir / "violin_prelim_violin.png"
        n_genes_plot = dataset_upload_dir / "scatter_prelim_n_genes.png"
        if violin_plot.is_file():
            try:
                shutil.move(violin_plot, dataset_final_dir / f"{dataset_id}.prelim_violin.png")
            except Exception as e:
                result['message'] = 'Error migrating preliminary violin plot: {}'.format(str(e))
                return result
        if n_genes_plot.is_file():
            try:
                shutil.move(n_genes_plot, dataset_final_dir / f"{dataset_id}.prelim_n_genes.png")
            except Exception as e:
                result['message'] = 'Error migrating preliminary n_genes plot: {}'.format(str(e))
                return result

    else:
        result['primary_analysis_migrated'] = 1

    # if we made it this far, all is well, so return success
    result['success'] = 1
    result['message'] = 'All steps completed successfully.'

    # now delete the entire upload directory
    shutil.rmtree(dataset_upload_dir)

    return result

if __name__ == '__main__':
    result = main()
    print(json.dumps(result))
