#!/opt/bin/python3

"""
Used by the expression uploader's "Review column types" step.

Applies the user's chosen dtype (categorical vs continuous) for each
flagged obs column directly onto the staged H5AD/Zarr file, then marks the
dataset as reviewed in metadata.json. This is the one point in the review
step where the dataset file is opened -- and even here, only obs is
touched (X, images, points, shapes, etc. are left alone).

Reuses the same adapters and read/mutate/write pattern already
established in gear.analysis / gear.primary_analysis:
    - Non-spatial: H5adAdapter(...).get_adata(backed=True), mutate .obs,
      adata.write() (writes back to the same path the backed object was
      opened from).
    - Spatial: ZarrAdapter(...).get_adata() returns sdata.tables["table"]
      (a plain AnnData); mutate .obs, then write just that table's zarr
      group back with adata.write_zarr(), leaving images/points/shapes
      untouched.

Expects a POST with:
    session_id
    share_uid
    choices  -- JSON string, e.g. '{"replicate": "categorical", "slide_id": "categorical"}'

Updates the staged dataset file in:
    ../uploads/files/<session_id>/<share_uid>/
and: ../uploads/files/<session_id>/<share_uid>/metadata.json
"""

import cgi
import json
import sys
from pathlib import Path

gear_root = Path(__file__).resolve().parents[2]
lib_path = gear_root / 'lib'
sys.path.append(str(lib_path))
import geardb
from gear.analysis import H5adAdapter, ZarrAdapter
from gear.utils import apply_obs_dtype_choices

user_upload_file_base = '../uploads/files'


def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    share_uid = form.getvalue('share_uid')
    choices_raw = form.getvalue('choices')

    result = {'success': 0, 'message': ''}

    if not session_id or not share_uid or not choices_raw:
        result['message'] = 'session_id, share_uid, and choices are required.'
        print(json.dumps(result))
        return

    user = geardb.get_user_from_session_id(session_id)
    if not user:
        result['message'] = 'Invalid session_id.'
        print(json.dumps(result))
        return

    try:
        choices = json.loads(choices_raw)
    except json.JSONDecodeError as e:
        result['message'] = f'Invalid choices JSON: {e}'
        print(json.dumps(result))
        return

    dataset_dir = Path(user_upload_file_base) / session_id / share_uid
    metadata_file = dataset_dir / 'metadata.json'

    if not metadata_file.is_file():
        result['message'] = 'No metadata JSON file found for this dataset.'
        print(json.dumps(result))
        return

    with open(metadata_file, 'r') as f:
        metadata = json.load(f)

    dataset_format = metadata.get('dataset_format', '')

    try:
        if dataset_format == 'spatial':
            zarr_path = dataset_dir / f"{share_uid}.zarr"
            if not zarr_path.is_dir():
                raise FileNotFoundError(f"Spatial dataset file not found: {zarr_path}")

            adapter = ZarrAdapter(zarr_path)
            adata = adapter.get_adata()  # sdata.tables["table"]
            adata.obs = apply_obs_dtype_choices(adata.obs, choices)

            # Write back just this table's zarr group -- not the whole
            # SpatialData store, so images/points/shapes are untouched.
            table_path = zarr_path / "tables" / "table"
            adata.write_zarr(table_path)
        else:
            h5ad_path = dataset_dir / f"{share_uid}.h5ad"
            if not h5ad_path.is_file():
                raise FileNotFoundError(f"Dataset file not found: {h5ad_path}")

            adapter = H5adAdapter(h5ad_path)
            adata = adapter.get_adata(backed=True)
            adata.obs = apply_obs_dtype_choices(adata.obs, choices)
            try:
                adata.write()  # writes back to h5ad_path
            finally:
                adata.file.close()
    except Exception as e:
        result['message'] = f'Error applying dtype choices: {e}'
        print(f'ObsDtypeReview error: {e}', file=sys.stderr)
        print(json.dumps(result))
        return

    # Mark reviewed and clear the flagged-column list so this step doesn't
    # get re-prompted if the user revisits the uploader for this dataset.
    metadata['obs_dtype_reviewed'] = True
    metadata['questionable_obs_columns'] = {}
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=4)

    result['success'] = 1
    result['message'] = 'Column types updated.'
    print(json.dumps(result))


if __name__ == '__main__':
    main()
