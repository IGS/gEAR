#!/opt/bin/python3

"""
Used by the expression uploader's post-processing review step.

Applies the user's choices to the staged H5AD/Zarr file, then records them in metadata.json:
    - the dtype (categorical vs continuous) for each flagged obs column
    - the dataset's primary variable (#435): the obs column the dataset is mainly about, saved as
      adata.uns["gear_obs_columns"]["primary_variable"] for building default displays

Both are applied in one request so only one process writes the staged file. Only obs and uns are
touched (X, raw, images, points, shapes, etc. are left alone):
    - Non-spatial: obs and uns/gear_obs_columns are rewritten in place with h5py and
      anndata's write_elem (as gear.primary_analysis does). Writing a backed AnnData with
      adata.write() would delete the raw group.
    - Spatial: ZarrAdapter(...).get_adata() returns sdata.tables["table"] (a plain AnnData);
      obs and uns are updated, then just that table's zarr group is written back with
      adata.write_zarr(), leaving images/points/shapes untouched.

Expects a POST with:
    session_id
    share_uid
    choices           -- optional JSON string, e.g. '{"replicate": "categorical", "slide_id": "categorical"}'
    primary_variable  -- optional obs column name; an empty value removes the primary variable.
                         Must be a categorical column once the choices are applied.
At least one of choices and primary_variable must be given.

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
from gear.analysis import ZarrAdapter
from gear.ontology import OBS_COLUMN_MAP_KEY
from gear.utils.obs import apply_obs_dtype_choices, is_categorical_obs_column
from werkzeug.utils import secure_filename

user_upload_file_base = '../uploads/files'

PRIMARY_VARIABLE_KEY = 'primary_variable'


def primary_variable_error(obs, primary_variable):
    """Return why primary_variable can't be used as the primary variable, or None if it can."""
    if not primary_variable:
        return None
    if primary_variable not in obs.columns:
        return f'The primary variable "{primary_variable}" is not a column in the dataset.'
    if not is_categorical_obs_column(obs[primary_variable]):
        return f'The primary variable "{primary_variable}" must be a categorical column.'
    return None


def update_obs_column_map(obs_column_map, primary_variable):
    """Return the gear_obs_columns map with the primary variable set, or removed if it's empty."""
    obs_column_map = dict(obs_column_map or {})
    obs_column_map.pop(PRIMARY_VARIABLE_KEY, None)
    if primary_variable:
        obs_column_map[PRIMARY_VARIABLE_KEY] = primary_variable
    return obs_column_map


def apply_to_h5ad(h5ad_path, choices, primary_variable, set_primary_variable):
    """Apply the choices to an h5ad in place. Returns an error message, or None on success."""
    import h5py
    from anndata.io import read_elem, write_elem

    with h5py.File(h5ad_path, 'a') as h5:
        obs = apply_obs_dtype_choices(read_elem(h5['obs']), choices)
        error = primary_variable_error(obs, primary_variable)
        if error:
            return error

        if choices:
            write_elem(h5, 'obs', obs)
        if set_primary_variable:
            uns = h5.require_group('uns')
            existing = read_elem(uns[OBS_COLUMN_MAP_KEY]) if OBS_COLUMN_MAP_KEY in uns else {}
            write_elem(uns, OBS_COLUMN_MAP_KEY, update_obs_column_map(existing, primary_variable))
    return None


def apply_to_zarr(zarr_path, choices, primary_variable, set_primary_variable):
    """Apply the choices to a spatial upload's table. Returns an error message, or None on success."""
    adapter = ZarrAdapter(zarr_path)
    adata = adapter.get_adata()  # sdata.tables["table"]
    adata.obs = apply_obs_dtype_choices(adata.obs, choices)
    error = primary_variable_error(adata.obs, primary_variable)
    if error:
        return error

    if set_primary_variable:
        adata.uns[OBS_COLUMN_MAP_KEY] = update_obs_column_map(adata.uns.get(OBS_COLUMN_MAP_KEY), primary_variable)

    # Write back just this table's zarr group -- not the whole
    # SpatialData store, so images/points/shapes are untouched.
    adata.write_zarr(zarr_path / "tables" / "table")
    return None


def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = secure_filename(form.getfirst('session_id', ''))
    share_uid = secure_filename(form.getfirst('share_uid', ''))
    choices_raw = form.getfirst('choices')
    # Absent leaves the primary variable as it is; empty removes it
    set_primary_variable = 'primary_variable' in form
    primary_variable = (form.getfirst('primary_variable') or '').strip() or None

    result = {'success': 0, 'message': ''}

    if not session_id or not share_uid or (not choices_raw and not set_primary_variable):
        result['message'] = 'session_id, share_uid, and choices or primary_variable are required.'
        print(json.dumps(result))
        return

    user = geardb.get_user_from_session_id(session_id)
    if not user:
        result['message'] = 'Invalid session_id.'
        print(json.dumps(result))
        return

    try:
        choices = json.loads(choices_raw) if choices_raw else {}
    except json.JSONDecodeError as e:
        result['message'] = f'Invalid choices JSON: {e}'
        print(json.dumps(result))
        return

    dataset_dir = (Path(user_upload_file_base) / session_id / share_uid).resolve()
    uploads_base = Path(user_upload_file_base).resolve()
    if not dataset_dir.is_relative_to(uploads_base):
        result['message'] = 'Invalid session_id or share_uid.'
        print(json.dumps(result))
        return

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
                raise FileNotFoundError(f"Spatial dataset file not found for share ID {share_uid}")
            error = apply_to_zarr(zarr_path, choices, primary_variable, set_primary_variable)
        else:
            h5ad_path = dataset_dir / f"{share_uid}.h5ad"
            if not h5ad_path.is_file():
                raise FileNotFoundError(f"Dataset file not found for share ID {share_uid}")
            error = apply_to_h5ad(h5ad_path, choices, primary_variable, set_primary_variable)
    except Exception as e:
        result['message'] = f'Error applying dtype choices: {e}'
        print(f'ObsDtypeReview error: {e}', file=sys.stderr)
        print(json.dumps(result))
        return

    if error:
        # Nothing was written
        result['message'] = error
        print(json.dumps(result))
        return

    # Mark reviewed and clear the flagged-column list so this step doesn't
    # get re-prompted if the user revisits the uploader for this dataset.
    if choices_raw:
        metadata['obs_dtype_reviewed'] = True
        metadata['questionable_obs_columns'] = {}
    if set_primary_variable:
        metadata['primary_variable'] = primary_variable
        metadata['primary_variable_reviewed'] = True
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=4)

    result['success'] = 1
    result['message'] = 'Column types updated.'
    print(json.dumps(result))


if __name__ == '__main__':
    main()
