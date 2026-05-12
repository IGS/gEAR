"""

rfuncs.py - Miscellaneous R-style functions called through rpy2
"""

import sys  # for print debugging
import traceback

import pandas as pd  # for dataframe manipulation

# rpy2.robjects calls rinterface.initr() under-the-hood when initialized to start the R session
import rpy2.robjects as ro

# If running locally, need to ensure that multiple concurrent R calls do not conflict
from rpy2.rinterface_lib import openrlib

# The number of R sessions appears to be limited to the number of threads apache allocates to the Flask API
# If this number of sessions exceeds number of threads, a RNotReady error will be thrown for each subsequent session
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr


class RError(Exception):
    """Error based on issues that would manifest in any particular R-language call."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)

def convert_r_df_to_r_matrix(df: pd.DataFrame):
    """
    Convert R-style dataframe to R-style matrix
    """

    r_matrix = ro.r["as.matrix"]
    return r_matrix(df, drop=False) # type: ignore

def convert_r_matrix_to_r_df(mtx) -> pd.DataFrame:
    """
    Convert R-style matrix to R-style dataframe
    """
    # mtx is a matrix of numbers with PCs in columns
    r_df = ro.r["as.data.frame"]
    return r_df(mtx) # type: ignore

def run_projectR_cmd(target_df: pd.DataFrame, loading_df: pd.DataFrame, algorithm: str, full_output:bool =False)-> list[pd.DataFrame]:
    """
    Convert input Pandas dataframes to R matrix.
    Pass the inputs into the projectR function written in R.
    Return Pandas dataframe of the projectR output
    """

    # Ensure multithreading if running locally -> https://rpy2.github.io/doc/v3.5.x/html/rinterface.html#multithreading
    with openrlib.rlock:

        # Convert from pandas dataframe to R data.frame
        # Seems any R call needs to be in a "conversion" context block
        # source -> https://stackoverflow.com/a/76532346, https://github.com/rpy2/rpy2/issues/1081, https://github.com/rpy2/rpy2/issues/975
        local_rules = localconverter(ro.default_converter + pandas2ri.converter)
        with local_rules:
            target_r_df = ro.conversion.py2rpy(target_df)
            loading_r_df = ro.conversion.py2rpy(loading_df)

            target_r_index = ro.conversion.py2rpy(target_df.index)
            loading_r_index = ro.conversion.py2rpy(loading_df.index)
        # Need a ruleset without pandas with auto-converts the R matrix to a numpy array

        with localconverter(ro.default_converter):
            try:
                # data.frame to data.matrix (projectR has no data.frame signature)
                target_r_matrix = convert_r_df_to_r_matrix(target_r_df)
                loading_r_matrix = convert_r_df_to_r_matrix(loading_r_df)
                # Assign Rownames to each matrix
                target_r_matrix.rownames = target_r_index
                loading_r_matrix.rownames = loading_r_index
            except Exception as e:
                # print stacktrace with line numbers
                traceback.print_exc(file=sys.stderr)
                raise RError("Error: Could not assign rownames to matrix.\tReason: {}".format(str(e)))

            # The NMF projectR method signature is based on the LinearEmbeddedMatrix class,
            # Which has a featureLoadings property. That matrix is loaded and the default
            # projectR signature is returned and used. So we can just pass the matrix as-is.
            # https://rdrr.io/bioc/SingleCellExperiment/man/LinearEmbeddingMatrix.html
            # Run project R command.  Get projectionPatterns matrix
            try:
                if algorithm == "nmf":
                    projectR = importr('projectR')
                    if full_output:
                        # R code: projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix)
                        projection_fit_r = projectR.projectR(data=target_r_matrix, loadings=loading_r_matrix, full=True)
                        # convert obj back to Python
                        projection_fit = ro.conversion.rpy2py(projection_fit_r)

                        # Both projection and pval are R-style matrices that need to be converted to R-style dataframes
                        projection_patterns_r_matrix = projection_fit[0]
                        pval_r_matrix = projection_fit[1]
                        projection_patterns_r_df = convert_r_matrix_to_r_df(projection_patterns_r_matrix)
                        pval_r_df = convert_r_matrix_to_r_df(pval_r_matrix)

                        with local_rules:
                            projection_patterns_df = ro.conversion.rpy2py(projection_patterns_r_df)
                            pval_df = ro.conversion.rpy2py(pval_r_df)
                            return [projection_patterns_df, pval_df]
                    else:
                        projection_patterns_r_matrix = projectR.projectR(data=target_r_matrix, loadings=loading_r_matrix, full=False)
                elif algorithm == "fixednmf":
                    sjd = importr('SJD')
                    loading_list = ro.ListVector({"genesig": loading_r_matrix})
                    projection = sjd.projectNMF(proj_dataset=target_r_matrix, proj_group=True, list_component=loading_list)
                    projection_patterns_r_matrix = projection.rx2("proj_score_list").rx2("genesig")
                else:
                    raise ValueError("Algorithm {} is not supported".format(algorithm))
            except ValueError as ve:
                # print stacktrace with line numbers
                traceback.print_exc(file=sys.stderr)
                raise
            except Exception as e:
                # print stacktrace with line numbers
                traceback.print_exc(file=sys.stderr)
                raise RError("Error: Could not run projectR command.\tReason: {}".format(str(e)))

            # matrix back to data.frame
            projection_patterns_r_df = convert_r_matrix_to_r_df(projection_patterns_r_matrix)

        with local_rules:
            # Convert from R data.frame to pandas dataframe
            projection_patterns_df = ro.conversion.rpy2py(projection_patterns_r_df)

            return [projection_patterns_df]


