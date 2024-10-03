#!/opt/bin/python3

# This script is used to update the share id (permalink)
# of either a layout, genecart, or dataset.

import cgi, json
import os
import sys
lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

from pathlib import Path
abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")
BY_DATASET_DIR = abs_path_www.joinpath("projections", "by_dataset")
BY_GENECART_DIR = abs_path_www.joinpath("projections", "by_genecart")


def main():
    cnx = geardb.Connection()
    print('Content-Type: application/json\n\n')

    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    share_id = form.getvalue('share_id')
    new_share_id = form.getvalue('new_share_id')
    scope = form.getvalue('scope') # 'dataset', 'layout', 'genecart'
    result = { 'error':"", 'success': 0 }

    user = geardb.get_user_from_session_id(session_id)
    if user is None:
        error = "Invalid session_id. User not found"
        result['error'] = error
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
        if gene_cart.gctype == "weighted":
            os.chdir(str(CARTS_BASE_DIR))

            for filename in os.listdir("."):
                if not share_id in filename:
                    continue
                if filename.endswith(".h5ad") or filename.endswith(".tab"):
                    # Replace old_id with new_id in the filename
                    new_filename = filename.replace(share_id, new_share_id)
                    #! This will not work if not owner or if permissions are not set correctly
                    os.rename(filename, new_filename)
                    # if the old file still exists, log it
                    if os.path.exists(filename):
                        print("Could not rename " + filename + " to " + new_filename + " for some reason... skipping", file=sys.stderr)

        # rename projections
        # the "by_dataset" directory only needs projections.json files updated
        os.chdir(str(BY_DATASET_DIR))
        for root, dirs, files in os.walk("."):
            for file in files:
                if not file.endswith("projections.json"):
                    continue
                filepath = os.path.join(root, file)
                with open(filepath, "r") as f:
                    lines = f.readlines()
                with open(filepath, "w") as f:
                    for line in lines:
                        f.write(line.replace(share_id, new_share_id))
                break

        # the "by_genecart" directory needs both directory names updated
        os.chdir(str(BY_GENECART_DIR))
        for root, dirs, files in os.walk("."):
            for dirname in dirs:
                if not share_id in dirname:
                    continue
                new_dir = dirname.replace(share_id, new_share_id)
                try:
                    os.rename(dirname, new_dir)
                except:
                    # If the new_dir already exists, we can't rename the directory
                    # The "by_genecart" directory is not exactly used, so not overly worried.
                    # But projections.json may need to be merged and files moved.
                    print("Could not rename " + dirname + " to " + new_dir + " because new_dir was not empty... skipping", file=sys.stderr)
                break

    else:
        error = "Invalid scope. Must be 'dataset', 'layout', or 'genecart'."
        result['error'] = error
        print(json.dumps(result))
        return

    # indicate success
    result['success'] = 1
    print(json.dumps(result))

if __name__ == '__main__':
    main()