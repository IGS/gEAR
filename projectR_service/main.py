import json, os
import pandas as pd
from flask import Flask, request

from rfuncs import run_projectR_cmd

app = Flask(__name__)

"""
@app.route("/")
def hello_world():
    name = os.environ.get("NAME", "World")
    return "Hello {}!".format(name)
"""

def decode_dataframe():
    pass

def encode_dataframe():
    pass

@app.route("/", methods=["POST"])
def index():

    target = request.form['target']
    loadings = request.form['loadings']
    is_pca = request.form['is_pca']

    target_df = pd.read_json(target)
    loading_df = pd.read_json(loadings)

    projection_patterns_df = run_projectR_cmd(target_df, loading_df, is_pca).transpose()

    return json.loads(projection_patterns_df.to_json())


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=int(os.environ.get("PORT", 8080)))