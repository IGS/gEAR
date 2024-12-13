import os, sys
import pandas as pd
from io import StringIO
from flask import Flask, abort, jsonify, request

cloud_logging = False
try:
    # Imports the Google Cloud client library
    from google.cloud import logging
    cloud_logging = True
except:
    pass


debug = os.environ.get('DEBUG', False)

app = Flask(__name__)

def write_entry(logger_name, severity, message):
    """Writes log entries to the given logger."""

    global cloud_logging

    if (not cloud_logging):
        print(message, file=sys.stderr)
        return

    logging_client = logging.Client()

    # This log can be found in the Cloud Logging console under 'Custom Logs'.
    logger = logging_client.logger(logger_name)

    # Simple text log with severity.
    logger.log_text(message, severity=severity)

### Each projection needs to return samples as rows, pattern weights as columns.

def do_binary_projection(target_df, loading_df):
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
        good_genes = target_df.loc[good_loading_genes].astype(bool).sum(axis=0).transpose()
        binary_target_df[pattern] = good_genes
    return binary_target_df


def do_pca_projection(target_df, loading_df):
    """Perform projection of PCA loadings."""
    tp_target_df = target_df.transpose()
    return tp_target_df.dot(loading_df)

@app.route("/status", methods=["GET"])
def status():
    return "OK"

@app.route("/", methods=["POST"])
def index():
    req_json = request.get_json()
    target = req_json['target']
    loadings = req_json['loadings']
    algorithm = req_json['algorithm']
    genecart_id = req_json["genecart_id"]
    dataset_id = req_json["dataset_id"]

    global cloud_logging
    try:
        write_entry("projectr", "DEBUG", "Testing cloud logging.")
    except Exception as e:
        cloud_logging = False
        write_entry("projectr", "DEBUG", "Failed to write to cloud logging: {}".format(str(e)))
        write_entry("projectr", "DEBUG", "Falling back to stderr logging.")

    write_entry("projectr", "INFO", "Dataset ID: {}".format(dataset_id))
    write_entry("projectr", "INFO", "Genecart ID: {}".format(genecart_id))


    # pd.read_json gives a FutureWarning, and suggest to wrap the json in StringIO.  Needed for pandas 2.x
    target = StringIO(target)
    loadings = StringIO(loadings)

    target_df = pd.read_json(target, orient="split")
    loading_df = pd.read_json(loadings, orient="split")

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

    # https://github.com/IGS/gEAR/issues/442#issuecomment-1317239909
    # Basically this is a stopgap until projectR has an option to remove
    # centering around zero for PCA loadings.  Chunking the data breaks
    # the output due to the centering around zero step.
    try:
        if algorithm == "pca":
            return jsonify(do_pca_projection(target_df,loading_df).to_json(orient="split"))
        elif algorithm == "binary":
            return jsonify(do_binary_projection(target_df, loading_df).to_json(orient="split"))
        elif algorithm == "2silca":
            pass
        elif algorithm in ["nmf", "fixednmf"]:
            from rfuncs import run_projectR_cmd
            projection_patterns_df = run_projectR_cmd(target_df, loading_df, algorithm).transpose()
        else:
            raise ValueError("Algorithm {} is not supported".format(algorithm))
    except Exception as e:
        description = str(e)
        write_entry("projectr", "ERROR", description)

        return abort(500, description=description)

    return jsonify(projection_patterns_df.to_json(orient="split"))

if __name__ == "__main__":
    app.run(debug=debug, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))