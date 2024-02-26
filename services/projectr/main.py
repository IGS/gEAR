import os, sys
import pandas as pd
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
    logging_client = logging.Client()

    # This log can be found in the Cloud Logging console under 'Custom Logs'.
    logger = logging_client.logger(logger_name)

    # Simple text log with severity.
    logger.log_text(message, severity=severity)

### Each projection needs to return samples as rows, pattern weights as columns.

def do_binary_projection(target_df, loading_df):
    """Perform projection based on the number of genes that were expressed in the cell or observation."""
    # Only applies with unweighted gene carts.
    tp_target_series = target_df.astype(bool).sum(axis=0).transpose()
    return pd.DataFrame(data=tp_target_series, columns=loading_df.columns, index=tp_target_series.index)

def do_pca_projection(target_df, loading_df):
    """Perform projection of PCA loadings."""
    tp_target_df = target_df.transpose()
    return tp_target_df.dot(loading_df)

@app.route("/", methods=["POST"])
def index():
    req_json = request.get_json()
    target = req_json['target']
    loadings = req_json['loadings']
    algorithm = req_json['algorithm']
    genecart_id = req_json["genecart_id"]
    dataset_id = req_json["dataset_id"]

    if cloud_logging:
        write_entry("projectr", "INFO", "Dataset ID: {}".format(dataset_id))
        write_entry("projectr", "INFO", "Genecart ID: {}".format(genecart_id))
    else:
        print("Dataset ID: {}".format(dataset_id), file=sys.stderr)
        print("Genecart ID: {}".format(genecart_id), file=sys.stderr)

    target_df = pd.read_json(target, orient="split")
    loading_df = pd.read_json(loadings, orient="split")

    if target_df.empty:
        description = "Target (dataset) dataframe is empty."
        if cloud_logging:
            write_entry("projectr", "ERROR", description)
        else:
            print(description, file=sys.stderr)
        return abort(500, description=description)

    if loading_df.empty:
        description = "Loading (pattern) dataframe is empty."
        if cloud_logging:
            write_entry("projectr", "ERROR", description)
        else:
            print(description, file=sys.stderr)
        return abort(500, description=description)

    if cloud_logging:
        write_entry("projectr", "INFO", "TARGET_DF SHAPE - {}".format(target_df.shape))
        write_entry("projectr", "INFO", "LOADING_DF SHAPE - {}".format(loading_df.shape))
    else:
        print("TARGET_DF SHAPE - {}".format(target_df.shape), file=sys.stderr)
        print("LOADING_DF SHAPE - {}".format(loading_df.shape), file=sys.stderr)

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
        if cloud_logging:
            write_entry("projectr", "ERROR", description)
        else:
            print(description, file=sys.stderr)
        return abort(500, description=description)

    return jsonify(projection_patterns_df.to_json(orient="split"))

if __name__ == "__main__":
    app.run(debug=debug, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))