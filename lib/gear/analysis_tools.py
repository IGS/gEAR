"""
analysis_tools.py - Which gEAR analysis tools support which dataset types (dtype).

Used by the availableAnalysisTools API (dataset explorer links) and by get_h5ad_dataset_list.cgi
(each tool's dataset list), so both agree on which datasets a tool can open. The lists keep legacy
dtypes (e.g. "sc-rna-seq") that older datasets still use.
"""

ANALYSIS_TOOL_DTYPES = {
    "dataset-curator": [
        "single-cell-rnaseq",
        "single-nucleus-rnaseq",
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
        "single-nucleus-rnaseq",
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
        "single-nucleus-rnaseq",
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
        "single-nucleus-rnaseq",
        "atac-seq",
        "spatial",
        "sc-rna-seq"
    ]
}
