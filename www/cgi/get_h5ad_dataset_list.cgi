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

    # allows for filtering datasets appropriate to specific pages
    for_page = form.getvalue("for_page")

    if session_id is not None:
        user = geardb.get_user_from_session_id(session_id)
    else:
        user = None

    public_collection = geardb.DatasetCollection()
    shared_with_user_collection = geardb.DatasetCollection()
    user_collection = geardb.DatasetCollection()

    result = {"user": user_collection, "public": public_collection, "shared_with_user": shared_with_user_collection}

    # public datasets are displayed no matter what
    if for_page is None:
        public_collection.get_public(has_h5ad=1)
    elif for_page == "analyze_dataset":
        public_collection.get_public(has_h5ad=1, types=["single-cell-rnaseq"])
    elif for_page == "compare_dataset":
        public_collection.get_public(
            has_h5ad=1,
            types=[
                "microarray",
                "bulk-rnaseq",
                "singlecell-h5ad",
                "single-cell-rnaseq",
                "svg-expression",
                "violin-standard",
            ],
        )
    elif for_page == "projection":
        public_collection.get_public(
            has_h5ad=1,
            types=[
                "microarray",
                "bulk-rnaseq",
                "singlecell-h5ad",
                "single-cell-rnaseq",
            ],
        )

    result["public"] = public_collection

    if user is not None:
        if for_page is None:
            user_collection = user.datasets(has_h5ad=1)
            shared_with_user_collection.get_shared_with_user(has_h5ad=1, user=user)
        elif for_page == "analyze_dataset":
            user_collection = user.datasets(has_h5ad=1, types=["single-cell-rnaseq"])
            shared_with_user_collection.get_shared_with_user(
                has_h5ad=1,
                user=user,
                types=[
                    "microarray",
                    "bulk-rnaseq",
                    "singlecell-h5ad",
                    "single-cell-rnaseq",
                    "violin-standard",
                ],
            )
        elif for_page == "compare_dataset":
            user_collection = user.datasets(
                has_h5ad=1,
                types=[
                    "microarray",
                    "bulk-rnaseq",
                    "singlecell-h5ad",
                    "single-cell-rnaseq",
                    "svg-expression",
                    "violin-standard",
                ],
            )
            shared_with_user_collection.get_shared_with_user(
                has_h5ad=1,
                user=user,
                types=[
                    "microarray",
                    "bulk-rnaseq",
                    "singlecell-h5ad",
                    "single-cell-rnaseq",
                    "svg-expression",
                    "violin-standard",
                ],
            )
        elif for_page == "projection":
            user_collection = user.datasets(
                has_h5ad=1,
                types=[
                    "microarray",
                    "bulk-rnaseq",
                    "singlecell-h5ad",
                    "single-cell-rnaseq",
                ],
            )
            shared_with_user_collection.get_shared_with_user(
                has_h5ad=1,
                user=user,
                types=[
                    "microarray",
                    "bulk-rnaseq",
                    "singlecell-h5ad",
                    "single-cell-rnaseq",
                ],
            )

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
