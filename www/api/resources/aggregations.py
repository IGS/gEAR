import typing

import geardb
from flask import request
from flask_restful import Resource, reqparse

from .common import get_adata_from_analysis, get_adata_shadow, get_spatial_adata

if typing.TYPE_CHECKING:
    from pandas import DataFrame

parser = reqparse.RequestParser(bundle_errors=True)
parser.add_argument("analysis_id", type=str, required=False)
parser.add_argument("filters", type=dict, required=False, default={})


class Aggregations(Resource):
    """
    Resource for handling aggregation requests on dataset observations.

    Methods
    -------
    post(dataset_id: str) -> dict
        Handles POST requests to aggregate categorical columns in the dataset's observations.
        Retrieves the AnnData object for the given dataset and analysis, applies optional filters,
        and returns the count of observations for each category in each categorical column.
        Categories with zero observations are included in the results.
        If the dataset file is not found, returns an error message.

    Parameters
    ----------
    dataset_id : str
        The unique identifier for the dataset to aggregate.

    Returns
    -------
    dict
        A dictionary containing:
            - "success": 1 if successful, -1 if file not found.
            - "aggregations": List of aggregation results for each categorical column.
            - "total_count": Total number of observations after filtering.
            - "message": Error message if file not found.
    """

    def post(self, dataset_id: str) -> dict:
        session_id = request.cookies.get("gear_session_id", "")
        args = parser.parse_args()
        analysis_id = args.get("analysis_id", None)
        filters = args.get("filters", {})  # key is column name, value is list of values

        ds = geardb.get_dataset_by_id(dataset_id)
        if not ds:
            return {
                "success": -1,
                'message': "No dataset found with that ID"
            }
        is_spatial = ds.dtype == "spatial"

        try:
            if is_spatial:
                adata = get_spatial_adata(analysis_id, dataset_id, session_id)
            elif filters:
                adata = get_adata_from_analysis(analysis_id, dataset_id, session_id)
            else:
                adata = get_adata_shadow(analysis_id, dataset_id, session_id)
        except FileNotFoundError:
            return {
                "success": -1,
                "aggregations": [],
                'message': "No dataset file found."
            }
        except Exception as e:
            return {
                "success": -1,
                "aggregations": [],
                'message': str(e)
            }

        obs: "DataFrame" = adata.obs  # type: ignore

        columns = obs.columns.tolist()

        if "replicate" in columns:
            columns.remove("replicate")

        columns = [col for col in columns if not col.endswith("_colors")]

        # Filter only categorical columns
        categorical_columns = list(
            filter(lambda x: obs[x].dtype.name == "category", columns)
        )

        # Get original categories for each categorical column
        # This is needed because we want to include categories with 0 observations
        orig_categories = {}
        for col in categorical_columns:
            orig_categories[col] = obs[col].cat.categories.tolist()

        # Filter currently selected values
        if filters:
            filter_query = ""
            for col, values in filters.items():
                if col in categorical_columns:
                    filter_query += f"(`{col}`.isin({values})) & "
            filter_query = filter_query[:-3]
            filter_idx = list(obs.query(filter_query).index)
            import numpy as np

            adata = adata[np.array(filter_idx), :]

        # Get number of observations for each value in each categorical column
        aggregations = []
        for col in categorical_columns:
            cat_aggregations = {"name": col, "count": obs[col].count(), "items": []}
            # Items has "name" and "count"
            value_dict = obs[col].astype("category").value_counts()

            # Add categories with 0 observations
            for cat in orig_categories[col]:
                if cat not in value_dict:
                    value_dict[cat] = 0

            # Rename "nan" to "Data not available"
            if "nan" in value_dict:
                value_dict["Data not available"] = value_dict.pop("nan")

            for value, count in value_dict.items():
                cat_aggregations["items"].append({"name": value, "count": count})
            aggregations.append(cat_aggregations)

        return {"success": 1, "aggregations": aggregations, "total_count": obs.shape[0]}
