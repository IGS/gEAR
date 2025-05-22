import os
import sys
from io import StringIO

import pandas as pd
from flask import Flask, abort, jsonify, request, Response

cloud_logging = False
try:
    # Imports the Google Cloud client library
    from google.cloud import logging

    cloud_logging = True
except Exception:
    pass


debug_str = os.environ.get("DEBUG", "False")
debug = False
if debug_str.lower() == "true":
    debug = True

app = Flask(__name__)


# create a 4-character random string
def random_string(length: int = 4) -> str:
    import random
    import string

    return "".join(random.choices(string.ascii_letters + string.digits, k=length))

identifier = random_string(4)

def write_entry(logger_name: str, severity: str, message: str) -> None:
    """Writes log entries to the given logger."""

    global cloud_logging

    message = "{} - {}".format(identifier, message)

    if not cloud_logging:
        print(message, file=sys.stderr)
        return

    logging_client = logging.Client()

    # This log can be found in the Cloud Logging console under 'Custom Logs'.
    logger = logging_client.logger(logger_name)

    # Simple text log with severity.
    logger.log_text(message, severity=severity)


### Each projection needs to return samples as rows, pattern weights as columns.


def do_binary_projection(
    target_df: pd.DataFrame, loading_df: pd.DataFrame
) -> pd.DataFrame:
    """Perform projection based on the number of genes that were expressed in the cell or observation."""
    # Only applies with unweighted gene carts, or weighted carts with binary values.

    # for each loading pattern, count the number of genes that are expressed in the target
    # and return the count as the pattern weight.
    binary_target_df = pd.DataFrame()
    for pattern in loading_df.columns:
        # Count the number of genes that are 1 in the loading_df
        good_loading_genes_mask = loading_df[pattern].astype(bool)
        good_loading_genes = loading_df.index[good_loading_genes_mask]

        # Count the number of those genes that are 1 (expressed) in the target_df.
        good_genes = (
            target_df.loc[good_loading_genes].astype(bool).sum(axis=0).transpose()
        )
        binary_target_df[pattern] = good_genes
    return binary_target_df


def do_pca_projection(
    target_df: pd.DataFrame, loading_df: pd.DataFrame
) -> pd.DataFrame:
    """Perform projection of PCA loadings."""
    tp_target_df = target_df.transpose()
    return tp_target_df.dot(loading_df)


@app.route("/status", methods=["GET"])
def status() -> str:
    """Health check endpoint."""
    return "OK"


@app.route("/", methods=["POST"])
def index() -> Response:
    req_json = request.get_json()
    target = req_json["target"]
    loadings = req_json["loadings"]
    algorithm = req_json["algorithm"]
    genecart_id = req_json["genecart_id"]
    dataset_id = req_json["dataset_id"]
    full_output = req_json.get("full_output", False)

    global cloud_logging
    try:
        write_entry("projectr", "DEBUG", "Testing cloud logging.")
    except Exception as e:
        cloud_logging = False
        write_entry(
            "projectr", "DEBUG", "Failed to write to cloud logging: {}".format(str(e))
        )
        write_entry("projectr", "DEBUG", "Falling back to stderr logging.")

    # Print the request payload to stderr
    write_entry("projectr", "DEBUG", "Genecart ID: {}".format(genecart_id))
    write_entry("projectr", "DEBUG", "Dataset ID: {}".format(dataset_id))
    write_entry("projectr", "DEBUG", "Algorithm: {}".format(algorithm))
    write_entry("projectr", "DEBUG", "Full output: {}".format(full_output))

    # pd.read_json gives a FutureWarning, and suggest to wrap the json in StringIO.  Needed for pandas 2.x
    target = StringIO(target)
    loadings = StringIO(loadings)

    target_df = pd.read_json(target, orient="split")
    loading_df = pd.read_json(loadings, orient="split")

    # fill NaN values with 0
    # This should have been done in the data processing step, but just in case.
    target_df = target_df.fillna(0)
    loading_df = loading_df.fillna(0)

    if target_df.empty:
        description = "Target (dataset) dataframe is empty."
        write_entry("projectr", "ERROR", description)
        return abort(500, description=description)

    if loading_df.empty:
        description = "Loading (pattern) dataframe is empty."
        write_entry("projectr", "ERROR", description)

        return abort(500, description=description)

    write_entry("projectr", "INFO", "TARGET_DF SHAPE - {}".format(target_df.shape))
    write_entry("projectr", "INFO", "LOADING_DF SHAPE - {}".format(loading_df.shape))

    response: dict["str", pd.DataFrame] = {
        "projection": pd.DataFrame(),
        "pval": pd.DataFrame(),
    }

    # https://github.com/IGS/gEAR/issues/442#issuecomment-1317239909
    # Basically this is a stopgap until projectR has an option to remove
    # centering around zero for PCA loadings.  Chunking the data breaks
    # the output due to the centering around zero step.
    try:
        if algorithm == "pca":
            response["projection"] = do_pca_projection(target_df, loading_df)
        elif algorithm == "binary":
            response["projection"] = do_binary_projection(target_df, loading_df)
        elif algorithm == "2silca":
            pass
        elif algorithm in ["nmf", "fixednmf"]:
            from rfuncs import run_projectR_cmd

            # R code: projectionFit <- list('projection'=projectionPatterns, 'pval'=pval.matrix)
            # "projection" is a DataFrame (index 0)
            # "pval" is a matrix (index 1)

            projection_patterns = run_projectR_cmd(
                target_df, loading_df, algorithm, full_output
            )

            projection_patterns_df = projection_patterns[0].transpose()
            response["projection"] = projection_patterns_df

            if full_output and algorithm == "nmf":
                projection_pval_df = projection_patterns[1].transpose()
                response["pval"] = projection_pval_df

        else:
            raise ValueError("Algorithm {} is not supported".format(algorithm))
    except Exception as e:
        description = str(e)
        write_entry("projectr", "ERROR", description)

        return abort(500, description=description)

    if response["projection"].empty:
        description = "Projection dataframe is empty."
        write_entry("projectr", "ERROR", description)

        return abort(500, description=description)

    write_entry("projectr", "INFO", "Projection dataframe shape: {}".format(response["projection"].shape))

    # Convert the response to JSON
    response_json = {
        "projection": response["projection"].to_json(orient="split"),
        "pval": response["pval"].to_json(orient="split"),
    }

    return jsonify(response_json)


if __name__ == "__main__":
    write_entry("projectr", "INFO", "Starting projectR service.")
    app.run(debug=debug, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))
