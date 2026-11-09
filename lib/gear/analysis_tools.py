"""
analysis_tools.py - Which gEAR analysis tools support which dataset types (dtype).

Used by the availableAnalysisTools API (dataset explorer links) and by the single-cell workbench's
dataset list, so both agree on which datasets a tool can open.
"""

ANALYSIS_TOOL_DTYPES = {
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
    # Single-cell data only; bulk datasets have too few observations for the workbench steps (#1361)
    "sc-workbench": [
        "single-cell-rnaseq",
        "atac-seq",
        "spatial",
        "sc-rna-seq"
    ]
}


def tool_supports_dtype(tool, dtype):
    """True if the named analysis tool can open a dataset of this dtype."""
    return dtype in ANALYSIS_TOOL_DTYPES.get(tool, [])
