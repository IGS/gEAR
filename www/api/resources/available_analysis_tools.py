"""
available_analysis_tools.py - Report which gEAR analysis tools support a dataset's type.

Serves /h5ad/<share_uid>/availableAnalysisTools in www/api/api.py.
"""

import geardb
from flask_restful import Resource
from gear.analysis_tools import ANALYSIS_TOOL_DTYPES


class AvailableAnalysisTools(Resource):
    """Available Analysis Tools

    Returns
    -------
    dict
        The available analysis tools for the specified dataset
    """
    def get(self, share_uid):
        """
        Return a dict of tool name to boolean availability based on the dataset's dtype.
        """
        dtype = geardb.get_dtype_by_share_id(share_uid)
        if not dtype:
            return {
                "success": -1,
                'message': "No dataset found with that share ID"
            }

        # return dict with [tool_name]: <boolean> indicating availability
        available_tools = {}
        for tool_name, dtypes in ANALYSIS_TOOL_DTYPES.items():
            available_tools[tool_name] = dtype in dtypes

        return {
            "success": 1,
            "share_id": share_uid,
            "available_analysis_tools": available_tools,
        }
