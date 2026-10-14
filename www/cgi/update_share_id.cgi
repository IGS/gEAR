#!/opt/bin/python3

# This script is used to update the share id (permalink)
# of either a layout, genecart, or dataset.

"""
update_share_id.cgi - Change the share ID (permalink) of a dataset, layout, or gene cart owned by the user.

Input: session_id (required), share_id, new_share_id, scope ('dataset', 'layout', 'genecart').
Output: JSON {success, error}. Gene carts also get cart files and projection paths renamed.
"""

import cgi
import json
import os
import sys
from pathlib import Path

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb
from werkzeug.utils import secure_filename

def is_safe_id(value: str) -> bool:
    """
    Custom share_ids can be user-chosen and vary in length, so we can't enforce a
    fixed format. Instead reject anything secure_filename() would alter -- this
    closes off path-injection via new_share_id, which is used unsanitized in
    filename.replace()/os.rename() below for gene carts.
    """
    return bool(value) and secure_filename(value) == value

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

# Longest share_id each scope's database column can hold (www/js/common.v2.js SHARE_ID_MAX_LENGTH)
SHARE_ID_MAX_LENGTH = {"dataset": 50, "genecart": 50, "layout": 24}
BY_DATASET_DIR = abs_path_www.joinpath("projections", "by_dataset")
BY_GENECART_DIR = abs_path_www.joinpath("projections", "by_genecart")

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    session_id = form.getfirst('session_id')
    share_id = form.getfirst('share_id')
    new_share_id = form.getfirst('new_share_id')
    scope = form.getfirst('scope') # 'dataset', 'layout', 'genecart'
    result = { 'error':"", 'success': 0 }

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        error = "Invalid session_id. User not found"
        result['error'] = error
        print(json.dumps(result))
        return

    if share_id is None:
        error = "Invalid share_id."
        result['error'] = error
        print(json.dumps(result))
        return

    if new_share_id is None:
        error = "Invalid new_share_id."
        result['error'] = error
        print(json.dumps(result))
        return

    if not is_safe_id(share_id) or not is_safe_id(new_share_id):
        error = ("This is not a valid permalink. Use only letters, numbers, hyphens, underscores and periods, "
                 "and don't start or end with an underscore or period.")
        result['error'] = error
        print(json.dumps(result))
        return

    # A longer value would fail to save (the column would reject it), so explain the limit instead
    max_length = SHARE_ID_MAX_LENGTH.get(scope)
    if max_length and len(new_share_id) > max_length:
        result['error'] = f"This is not a valid permalink. It can be at most {max_length} characters long."
        print(json.dumps(result))
        return

    # ? Should I add the validate queries to the geardb module?

    if scope == "dataset":
        # Verify that the share_id exists
        # Otherwise we have nothing to match against
        dataset_id = geardb.get_dataset_id_from_share_id(share_id)
        if dataset_id is None:
            error = "Invalid share_id."
            result['error'] = error
            print(json.dumps(result))
            return

        dataset = geardb.get_dataset_by_id(dataset_id)

        if dataset is None:
            error = "Dataset not found."
            result['error'] = error
            print(json.dumps(result))
            return

        # Verify this user owns the dataset
        if dataset.owner_id != user.id:
            error = "User does not own this dataset."
            result['error'] = error
            print(json.dumps(result))
            return

        # New share_id should not be present in the database
        dataset2_id = geardb.get_dataset_id_from_share_id(new_share_id)
        if dataset2_id is not None:
            error = "New share_id already exists."
            result['error'] = error
            print(json.dumps(result))
            return

        # Update the share_id
        dataset.save_change("share_id", new_share_id)

    elif scope == "layout":
        layout = geardb.get_layout_by_share_id(share_id)
        if layout is None:
            error = "Invalid share_id."
            result['error'] = error
            print(json.dumps(result))
            return

        if layout.user_id != user.id:
            error = "User does not own this layout."
            result['error'] = error
            print(json.dumps(result))
            return

        layout2 = geardb.get_layout_by_share_id(new_share_id)
        if layout2 is not None:
            error = "New share_id already exists."
            result['error'] = error
            print(json.dumps(result))
            return

        # Update the share_id
        layout.save_change("share_id", new_share_id)

    elif scope == "genecart":
        gene_cart = geardb.get_gene_cart_by_share_id(share_id)
        if gene_cart is None:
            error = "Invalid share_id."
            result['error'] = error
            print(json.dumps(result))
            return

        if gene_cart.user_id != user.id:
            error = "User does not own this genecart."
            result['error'] = error
            print(json.dumps(result))
            return

        gene_cart2 = geardb.get_gene_cart_by_share_id(new_share_id)
        if gene_cart2 is not None:
            error = "New share_id already exists."
            result['error'] = error
            print(json.dumps(result))
            return

        gene_cart.save_change("share_id", new_share_id)

        # Gene carts are special cases.
        # If this is a weighted gene cart, we need to update the share_id used
        # in the filename of www/carts<share_id>.(h5ad|tab)
        # Any projection.json files that reference this genecart will also need to be updated.

        # rename carts
        if gene_cart.gctype == "weighted-list":
            os.chdir(str(CARTS_BASE_DIR))

            for ext in (".h5ad", ".tab"):
                # Only this cart's own files: another cart's name may contain share_id
                #  (e.g. cart.<share_id>99.tab), and must not be renamed with it
                filename = f"cart.{share_id}{ext}"
                if not os.path.exists(filename):
                    continue
                new_filename = f"cart.{new_share_id}{ext}"

                # Defense in depth: new_share_id was already validated by is_safe_id(),
                # but confirm the resulting path still resolves inside CARTS_BASE_DIR
                # before renaming.
                resolved_new_filename = (CARTS_BASE_DIR / new_filename).resolve()
                if not resolved_new_filename.is_relative_to(CARTS_BASE_DIR.resolve()):
                    print("Refusing to rename " + filename + " to " + new_filename + ": resolves outside carts directory", file=sys.stderr)
                    continue

                #! This will not work if not owner or if permissions are not set correctly
                os.rename(filename, resolved_new_filename)
                # if the old file still exists, log it
                if os.path.exists(filename):
                    print("Could not rename " + filename + " to " + new_filename + " for some reason... skipping", file=sys.stderr)

        # rename projections
        rename_projection_keys(BY_DATASET_DIR, share_id, new_share_id)
        rename_genecart_projection_dir(BY_GENECART_DIR, share_id, new_share_id)

    else:
        error = "Invalid scope. Must be 'dataset', 'layout', or 'genecart'."
        result['error'] = error
        print(json.dumps(result))
        return

    # indicate success
    result['success'] = 1
    print(json.dumps(result))

def rename_projection_keys(by_dataset_dir: Path, share_id: str, new_share_id: str) -> None:
    """
    Rename the gene cart's key in each dataset's projections.json. Only the key is changed: other
    gene carts, datasets and projection IDs in the file may contain share_id as a substring.
    """
    if not by_dataset_dir.is_dir():
        return

    old_keys = (share_id, "cart.{}".format(share_id))   # "cart.<id>" is the legacy key
    for json_file in by_dataset_dir.glob("*/projections.json"):
        try:
            with open(json_file) as fh:
                projections = json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            print("Could not read {}: {}... skipping".format(json_file, e), file=sys.stderr)
            continue

        if not isinstance(projections, dict) or not any(key in projections for key in old_keys):
            continue

        for key in old_keys:
            if key in projections:
                projections.setdefault(new_share_id, []).extend(projections.pop(key))

        with open(json_file, "w") as fh:
            json.dump(projections, fh, ensure_ascii=False, indent=4)

def rename_genecart_projection_dir(by_genecart_dir: Path, share_id: str, new_share_id: str) -> None:
    """Rename the gene cart's own directory under projections/by_genecart, if it has one."""
    old_dir = by_genecart_dir / share_id
    if not old_dir.is_dir():
        return

    new_dir = by_genecart_dir / new_share_id
    # Defense in depth: new_share_id was already validated by is_safe_id()
    if not new_dir.resolve().is_relative_to(by_genecart_dir.resolve()):
        print("Refusing to rename " + str(old_dir) + " to " + str(new_dir) + ": resolves outside genecart projections directory", file=sys.stderr)
        return
    if new_dir.exists():
        # The "by_genecart" directory is not exactly used, so not overly worried.
        # But projections.json may need to be merged and files moved.
        print("Could not rename " + str(old_dir) + " to " + str(new_dir) + " because it already exists... skipping", file=sys.stderr)
        return

    old_dir.rename(new_dir)

if __name__ == '__main__':
    main()