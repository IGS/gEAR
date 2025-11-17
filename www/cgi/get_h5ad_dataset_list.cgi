#!/opt/bin/python3

"""
Used by sc_workbench.html, this script gets a list of the H5AD datasets the user can view.
Returns two sets 'public' and 'user' (if a user can be pulled from the session)
"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join("..", "..", "lib"))
sys.path.append(lib_path)
import geardb


def main():
    form = cgi.FieldStorage()
    session_id = form.getvalue("session_id")
    share_id = form.getvalue("share_id", None)
    if share_id == "null":
        share_id = None

    if session_id is not None:
        user = geardb.get_user_from_session_id(session_id)
    else:
        user = None

    public_collection = geardb.DatasetCollection()
    shared_with_user_collection = geardb.DatasetCollection()
    user_collection = geardb.DatasetCollection()

    # Spatial datasets have no h5ad but have the "spatial" dtype.
    # Downstream, we can convert SpatialData objects to AnnData so we should include these.
    spatial_public_collection = geardb.DatasetCollection()
    spatial_shared_with_user_collection = geardb.DatasetCollection()
    spatial_user_collection = geardb.DatasetCollection()

    result = {"user": user_collection, "public": public_collection, "shared_with_user": shared_with_user_collection}

    # public datasets are displayed no matter what
    public_collection.get_public(has_h5ad=1)
    spatial_public_collection.get_public(has_h5ad=0, types=["spatial"])
    public_collection.datasets.extend(spatial_public_collection.datasets)

    result["public"] = public_collection

    if user is not None:
        # NOTE: This previous used user.datasets() but that does not clear the "_datasets" property if called multiple times.
        user_collection.get_owned_by_user(has_h5ad=1, user=user)
        spatial_user_collection.get_owned_by_user(has_h5ad=0, user=user, types=["spatial"])
        user_collection.datasets.extend(spatial_user_collection.datasets)

        shared_with_user_collection.get_shared_with_user(has_h5ad=1, user=user)
        spatial_shared_with_user_collection.get_shared_with_user(has_h5ad=0, types=["spatial"], user=user)
        shared_with_user_collection.datasets.extend(spatial_shared_with_user_collection.datasets)


    result["user"] = user_collection
    result["shared_with_user"] = shared_with_user_collection

    def dataset_in_collections(dataset_id, collections):
        return any(
            ds.id == dataset_id
            for colname in collections
            for ds in result[colname].datasets
        )

    # If a dataset_id was directly passed, check if it was included in any of the
    # collections. If not, manually add it.
    if share_id:
        included_ds = geardb.get_dataset_by_share_id(share_id=share_id)
        if not included_ds:
            print("Content-Type: application/json\n\n")
            print(json.dumps(result))
            return

        included_ds_id = included_ds.id
        collections = ["user", "public", "shared_with_user"]

        if not dataset_in_collections(included_ds_id, collections):
            result["shared_with_user"].datasets.append(included_ds)

    print("Content-Type: application/json\n\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
