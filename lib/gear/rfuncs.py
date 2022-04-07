"""

rfuncs.py - Miscellaneous R-style functions called through rpy2
"""

import sys  # for print debugging

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import ListVector, StrVector

def convert_df_to_r_matrix(df):
    """
    Convert pandas dataframe to R-style matrix
    """
    r_matrix = robjects.r["as.matrix"]
    return r_matrix(df)

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

    projectR = importr('projectR')

    # Convert from pandas dataframe to R data.frame
    with localconverter(robjects.default_converter + pandas2ri.converter):
       target_r_df = robjects.conversion.py2rpy(target_df)
       loading_r_df = robjects.conversion.py2rpy(loading_df)

    # data.frame to matrix (projectR has no data.frame signature)
    target_r_matrix = convert_df_to_r_matrix(target_r_df)
    loading_r_matrix = convert_df_to_r_matrix(loading_r_df)

    # Assign Rownames to each matrix
    target_r_matrix.rownames = StrVector(target_df.index)
    loading_r_matrix.rownames = StrVector(loading_df.index)

    loading_r_object = loading_r_matrix

    if is_pca:
        loading_r_object = convert_r_matrix_to_r_prcomp(loading_r_matrix)

    # Run project R command.  Get projectionPatterns matrix
    projection_patterns_r_matrix = projectR.projectR(data=target_r_matrix, loadings=loading_r_object, full=False)

    # matrix back to data.frame
    r_df = robjects.r["as.data.frame"]
    projection_patterns_r_df = r_df(projection_patterns_r_matrix)

    # Convert from R data.frame to pandas dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
       projection_patterns_df = robjects.conversion.rpy2py(projection_patterns_r_df)

    return projection_patterns_df


