import json, os, sys
import pandas as pd
from flask import Flask, request

from rfuncs import run_projectR_cmd, RError

app = Flask(__name__)

@app.route("/", methods=["POST"])
def index():
    req_json = request.get_json()
    target = req_json['target']
    loadings = req_json['loadings']
    is_pca = req_json['is_pca']
    genecart_id = req_json["genecart_id"]
    dataset_id = req_json["dataset_id"]

    # TODO: change print msgs to be structured log messages
    # https://cloud.google.com/run/docs/samples/cloudrun-manual-logging
    print("Dataset ID: {}".format(dataset_id), file=sys.stderr)
    print("Genecart ID: {}".format(genecart_id), file=sys.stderr)

    target_df = pd.read_json(target, orient="split")
    loading_df = pd.read_json(loadings, orient="split")


    if target_df.empty:
        raise RError("Target (dataset) dataframe is empty.")

    if loading_df.empty:
        raise RError("Loading (pattern) dataframe is empty.")

    projection_patterns_df = run_projectR_cmd(target_df, loading_df, is_pca).transpose()

    return json.loads(projection_patterns_df.to_json())


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))