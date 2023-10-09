from flask import request
from flask_restful import Resource
import os
import geardb


class Aggregations(Resource):
    """Resource for retrieving observation aggregations for a dataset and applied categorial observations filters

    Parameters
    ----------
    dataset_id: str
        Dataset ID
    session_id: str
        Session ID
    analysis_id: str
        Analysis ID
    filters: dict
        Filters applied to dataset. Key is column name, value is list of values

    Returns
    -------
    list of dicts
        * name: categorical column name
        * count: number of observations
        * items: list of dicts
            * name: categorical value
            * count: number of observations


    """

    def post(self, dataset_id):
        req = request.get_json()
        dataset_id = req.get('dataset_id')
        session_id = req.get('session_id')
        analysis_id = req.get('analysis_id')
        filters = req.get('filters')    # key is column name, value is list of values


        ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
        h5_path = ds.get_file_path()

        # Import here so that Flask-RESTful does not import it with every API call.
        import scanpy as sc

        # Have a public dataset or user_saved dataset
        if analysis_id:
            user = geardb.get_user_from_session_id(session_id)

            ana = geardb.Analysis(id=analysis_id, dataset_id=dataset_id, session_id=session_id, user_id=user.id)
            ana.discover_type()

        else:
            # Dataset is primary type
            if not os.path.exists(h5_path):
                # Let's not fail if the file isn't there
                return {
                    "success": -1,
                    'message': "No h5 file found for this dataset"
                }
            adata = sc.read_h5ad(h5_path)

        columns = adata.obs.columns.tolist()

        if "replicate" in columns:
            columns.remove('replicate')

        columns = [col for col in columns if not col.endswith('_colors')]

        # Filter only categorical columns
        categorical_columns = list(filter(lambda x: adata.obs[x].dtype.name == 'category', columns))

        # Get original categories for each categorical column
        # This is needed because we want to include categories with 0 observations
        orig_categories = {}
        for col in categorical_columns:
            orig_categories[col] = adata.obs[col].cat.categories.tolist()

        # Filter currently selected values
        if filters:
            filter_query = ""
            for col, values in filters.items():
                if col in categorical_columns:
                    filter_query += f"(`{col}`.isin({values})) & "
            filter_query = filter_query[:-3]
            adata = adata[adata.obs.query(filter_query).index]

        # Get number of observations for each value in each categorical column
        aggregations = []
        for col in categorical_columns:
            cat_aggregations = {"name": col, "count": adata.obs[col].count(),"items": []}
            # Items has "name" and "count"
            value_dict = adata.obs[col].astype("category").value_counts()

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

        return {
            "success": 1,
            "aggregations": aggregations,
            "total_count": adata.obs.shape[0]
        }


