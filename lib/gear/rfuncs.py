"""

rfuncs.py - Miscellaneous R-style functions called through rpy2
"""

import sys  # for print debugging
import traceback    # for debugging
import gc   # garbage collection

import rpy2.rinterface as ri    # Since the low-level interface does not auto-intialize R, we can add globally.
from rpy2.rinterface_lib import openrlib

from time import sleep

class RError(Exception):
    """Error based on issues that would manifest in any particular R-language call."""
    def __init__(self, message="") -> None:
        self.message = message
        super().__init__(self.message)

def convert_r_df_to_r_matrix(df):
    """
    Convert pandas dataframe to R-style matrix
    """
    r_matrix = ri.baseenv["as.matrix"]
    return r_matrix(df)

def convert_r_matrix_to_r_df(mtx):
    """
    Convert R-style matrix to R-style dataframe
    """
    # mtx is a matrix of numbers with PCs in columns
    r_df = ri.baseenv["as.data.frame"]
    return r_df(mtx)

def convert_r_matrix_to_r_prcomp(mtx):
    """
    Convert R-style matrix to R-style prcomp object
    """
    # mtx is a matrix of numbers with PCs in columns
    from rpy2.robjects import ListVector
    prcomp_obj = ListVector({"rotation": mtx})
    prcomp_obj.rclass = "prcomp"    # Convert to prcomp-class object
    return prcomp_obj

def run_projectR_cmd(target_df, loading_df, is_pca=False):
    """
    Convert input Pandas dataframes to R matrix.
    Pass the inputs into the projectR function written in R.
    Return Pandas dataframe of the projectR output
    """

    # Import the low-level R-interface first and perform a manual initialization of the R session
    # rpy2.robjects also calls initr() under-the-hood when imported but subsequent calls are ignored
    # The number of R sessions appears to be limited to the number of threads apache allocates to the Flask API
    # If this number of sessions exceeds number of threads, a RNotReady error will be thrown for each subsequent session
    #try:
    ri.initr_simple()
    #except Exception as e:
    #    print("ERROR: {}".format(str(e)), file=sys.stderr)
    #    ri.endr(1)  # Exit with fatal state
    #    raise RError("Could not initialize R session to run projectR.")
    #sleep(5)    # Give enough time for the R session to start

    # R does not play nice with multithreading so a lock is necessary to prevent interruptions
    with openrlib.rlock:
        if target_df.empty:
            #ri.endr(1)  # Exit with fatal state
            raise RError("Target (dataset) dataframe is empty.")

        if loading_df.empty:
            #ri.endr(1)  # Exit with fatal state
            raise RError("Loading (pattern) dataframe is empty.")

        # NOTE: Importing robjects inside of function so the Flask-RESTful API does not initialize R at the beginning of every API call
        #import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri, default_converter
        from rpy2.robjects.packages import importr
        from rpy2.robjects.conversion import localconverter, py2rpy, rpy2py

        #from rpy2.robjects.conversion import localconverter, py2rpy, rpy2py

        # Convert from pandas dataframe to R data.frame
        with localconverter(default_converter + pandas2ri.converter):
            target_r_df = py2rpy(target_df)
            loading_r_df = py2rpy(loading_df)

        # data.frame to matrix (projectR has no data.frame signature)
        target_r_matrix = convert_r_df_to_r_matrix(target_r_df)
        loading_r_matrix = convert_r_df_to_r_matrix(loading_r_df)

        # Assign Rownames to each matrix
        target_r_matrix.rownames = ri.ListSexpVector(target_df.index)    # low-level equivalent to ro.StrVector
        loading_r_matrix.rownames = ri.ListSexpVector(loading_df.index)

        # Modify the R-style matrix to be a prcomp object if necessary
        loading_r_object = loading_r_matrix

        if is_pca:
            loading_r_object = convert_r_matrix_to_r_prcomp(loading_r_matrix)
            if not loading_r_object:
                ri.endr(1)  # Exit with fatal state
                raise RError("Could not convert loading matrix to R prcomp object.")

        # Run project R command.  Get projectionPatterns matrix
        try:
            projectR = importr('projectR')
            projection_patterns_r_matrix = projectR.projectR(data=target_r_matrix, loadings=loading_r_object, full=False)
        except Exception as e:
            print("ERROR: {}".format(str(e)), file=sys.stderr)
            #print(traceback.print_exc(), file=sys.stderr)
            gc.collect()    # Reduce chances of memory issues
            ri.endr(1)  # Exit with fatal state
            raise RError("Error: Could not run projectR command.")

        # matrix back to data.frame
        projection_patterns_r_df = convert_r_matrix_to_r_df(projection_patterns_r_matrix)

        # Convert from R data.frame to pandas dataframe
        with localconverter(default_converter + pandas2ri.converter):
            projection_patterns_df = rpy2py(projection_patterns_r_df)

        gc.collect()
        ri.endr(0)
        return projection_patterns_df


