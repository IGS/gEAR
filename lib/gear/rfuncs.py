"""

rfuncs.py - Miscellaneous R-style functions called through rpy2
"""

import sys  # for print debugging
import traceback    # for debugging

import gc

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import ListVector, StrVector
from rpy2.rinterface_lib import openrlib

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

def convert_r_matrix_to_r_prcomp(mtx):
    """
    Convert R-style matrix to R-style prcomp object
    """
    # mtx is a matrix of numbers with PCs in columns
    prcomp_obj = ListVector({"rotation": mtx})
    prcomp_obj.rclass = "prcomp"    # Convert to prcomp-class object
    return prcomp_obj

def run_projectR_cmd(target_df, loading_df, is_pca=False):
    """
    Convert input Pandas dataframes to R matrix.
    Pass the inputs into the projectR function written in R.
    Return Pandas dataframe of the projectR output
    """

    # R does not play nice with multithreading so a lock is necessary to prevent interruptions
    try:
        projectR = importr('projectR')
    except Exception as e:
        print("ERROR: {}".format(str(e)), file=sys.stderr)
        print(traceback.print_exc(), file=sys.stderr)
        raise RError("Could not import projectR package.")

    # Convert from pandas dataframe to R data.frame
    with localconverter(ro.default_converter + pandas2ri.converter):
        target_r_df = ro.conversion.py2rpy(target_df)
        loading_r_df = ro.conversion.py2rpy(loading_df)

    # data.frame to matrix (projectR has no data.frame signature)
    target_r_matrix = convert_r_df_to_r_matrix(target_r_df)
    loading_r_matrix = convert_r_df_to_r_matrix(loading_r_df)

    # Assign Rownames to each matrix
    target_r_matrix.rownames = StrVector(target_df.index)
    loading_r_matrix.rownames = StrVector(loading_df.index)

    # Modify the R-style matrix to be a prcomp object if necessary
    loading_r_object = loading_r_matrix
    if is_pca:
        loading_r_object = convert_r_matrix_to_r_prcomp(loading_r_matrix)

    with openrlib.rlock:
        # Run project R command.  Get projectionPatterns matrix
        try:
            projection_patterns_r_matrix = projectR.projectR(data=target_r_matrix, loadings=loading_r_object, full=False)
        except Exception as e:
            gc.collect()    # collect R garbage
            print("ERROR: {}".format(str(e)), file=sys.stderr)
            print(traceback.print_exc(), file=sys.stderr)
            raise RError("Error: Could not run projectR command.")

    # matrix back to data.frame
    projection_patterns_r_df = convert_r_matrix_to_r_df(projection_patterns_r_matrix)

    # Convert from R data.frame to pandas dataframe
    with localconverter(ro.default_converter + pandas2ri.converter):
        projection_patterns_df = ro.conversion.rpy2py(projection_patterns_r_df)

    gc.collect()    # collect R garbage
    return projection_patterns_df


