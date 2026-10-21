"""Fake gear.metadata: Metadata.save_to_mysql logs to the fake database log instead of writing."""

import json

import geardb  # the fake geardb installed by run_cgi.py


class Metadata:
    def __init__(self, file_path=None, **kwargs):
        self.file_path = file_path
        with open(file_path) as fh:
            self.fields = json.load(fh)

    def make_spatial_h5ad_adjustment(self, dataset_format):
        geardb._log("metadata.make_spatial_h5ad_adjustment", dataset_format=dataset_format)

    def save_to_mysql(self, status=None, is_public=0):
        geardb._log("metadata.save_to_mysql", dataset_uid=self.fields.get("dataset_uid"),
                    status=status, is_public=is_public)


def get_value_from_df(df, key):
    return None
