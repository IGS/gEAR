import json, os, sys
import pandas as pd
from flask import Flask, request

from rfuncs import run_projectR_cmd

app = Flask(__name__)

@app.route("/", methods=["POST"])
def index():
    req_json = request.get_json()
    target = req_json['target']
    loadings = req_json['loadings']
    is_pca = req_json['is_pca']

    target_df = pd.read_json(target)
    loading_df = pd.read_json(loadings)

    projection_patterns_df = run_projectR_cmd(target_df, loading_df, is_pca).transpose()

    return json.loads(projection_patterns_df.to_json())


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))