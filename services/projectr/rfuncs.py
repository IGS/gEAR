"""

rfuncs.py - Miscellaneous R-style functions called through rpy2
"""

import sys  # for print debugging

# rpy2.robjects calls rinterface.initr() under-the-hood when initialized to start the R session
import rpy2.robjects as ro
# The number of R sessions appears to be limited to the number of threads apache allocates to the Flask API
# If this number of sessions exceeds number of threads, a RNotReady error will be thrown for each subsequent session

from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import ListVector, StrVector

class RError(Exception):
    """Error based on issues that would manifest in any particular R-language call."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)

def convert_r_df_to_r_matrix(df):
    """
    Convert pandas dataframe to R-style matrix
    """
    r_matrix = ro.r["as.matrix"]
    return r_matrix(df)

def convert_r_matrix_to_r_df(mtx):
    """
    Convert R-style matrix to R-style dataframe
    """
    # mtx is a matrix of numbers with PCs in columns
    r_df = ro.r["as.data.frame"]
    return r_df(mtx)

def run_projectR_cmd(target_df, loading_df):
    """
    Convert input Pandas dataframes to R matrix.
    Pass the inputs into the projectR function written in R.
    Return Pandas dataframe of the projectR output
    """

    # Convert from pandas dataframe to R data.frame
    with localconverter(ro.default_converter + pandas2ri.converter):
        target_r_df = ro.conversion.py2rpy(target_df)
        loading_r_df = ro.conversion.py2rpy(loading_df)

    # data.frame to matrix (projectR has no data.frame signature)
    target_r_matrix = convert_r_df_to_r_matrix(target_r_df)
    loading_r_matrix = convert_r_df_to_r_matrix(loading_r_df)

    # Assign Rownames to each matrix
    # I don't know why but using ro.StrVector makes rpy2py fail where the output df is an incompatible class
    # Guessing that there are some non-strings mixed into the indexes
    target_r_matrix.rownames = StrVector(target_df.index)
    loading_r_matrix.rownames = StrVector(loading_df.index)

    # The NMF projectR method signature is based on the LinearEmbeddedMatrix class,
    # Which has a featureLoadings property. That matrix is loaded and the default
    # projectR signature is returned and used. So we can just pass the matrix as-is.
    # https://rdrr.io/bioc/SingleCellExperiment/man/LinearEmbeddingMatrix.html

    # Run project R command.  Get projectionPatterns matrix
    try:
        projectR = importr('projectR')
        projection_patterns_r_matrix = projectR.projectR(data=target_r_matrix, loadings=loading_r_matrix, full=False)
    except Exception as e:
        raise RError("Error: Could not run projectR command.\tReason: {}".format(str(e)))

    # matrix back to data.frame
    projection_patterns_r_df = convert_r_matrix_to_r_df(projection_patterns_r_matrix)

    # Convert from R data.frame to pandas dataframe
    with localconverter(ro.default_converter + pandas2ri.converter):
        projection_patterns_df = ro.conversion.rpy2py(projection_patterns_r_df)

    return projection_patterns_df


