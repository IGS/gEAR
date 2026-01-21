import geardb
from flask_restful import Resource

tools = {
    "dataset-curator": [
        "single-cell-rnaseq",
        "bulk-rnaseq",
        "bargraph-standard",
        "microarray",
        "svg-expression",
        "atac-seq",
        "violin-standard",
        "spatial",
        "sc-rna-seq",
        "linegraph-standard"
    ],
    "multigene-viewer": [
        "single-cell-rnaseq",
        "bulk-rnaseq",
        "bargraph-standard",
        "microarray",
        "atac-seq",
        "violin-standard",
        "spatial",
        "sc-rna-seq",
        "linegraph-standard"
    ],
    "compare-tool": [
        "single-cell-rnaseq",
        "bulk-rnaseq",
        "bargraph-standard",
        "microarray",
        "atac-seq",
        "violin-standard",
        "spatial",
        "sc-rna-seq",
        "linegraph-standard"
    ],
    "sc-workbench": [
        "single-cell-rnaseq",
        "atac-seq",
        "spatial",
        "sc-rna-seq"
    ]
}

class AvailableAnalysisTools(Resource):
    """Available Analysis Tools

    Returns
    -------
    dict
        The available analysis tools for the specified dataset
    """
    def get(self, share_uid):
        dtype = geardb.get_dtype_by_share_id(share_uid)
        if not dtype:
            return {
                "success": -1,
                'message': "No dataset found with that share ID"
            }

        # return dict with [tool_name]: <boolean> indicating availability
        available_tools = {}
        for tool_name, dtypes in tools.items():
            available_tools[tool_name] = dtype in dtypes

        return {
            "success": 1,
            "share_id": share_uid,
            "available_analysis_tools": available_tools,
        }
