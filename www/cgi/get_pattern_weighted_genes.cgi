#!/opt/bin/python3

"""
For a given weighted genecart, return all genes and their weights (sorted descending) for the specified pattern

Input: source_id (gene cart share ID), pattern_id (weight column name).
Output: JSON list of {gene, weight}; on error, HTTP 400/404 with JSON {success: 0, error}.
"""

import cgi
import json
import pandas as pd
from pathlib import Path

from werkzeug.utils import secure_filename

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

def print_error(status, message):
    """Print a JSON error response with the given HTTP status."""
    print(f"Status: {status}")
    print('Content-Type: application/json\n\n')
    print(json.dumps({'success': 0, 'error': message}))

def main():
    form = cgi.FieldStorage()
    source_id = form.getfirst('source_id')  # Root of the file name (minus extension)
    pattern_id = form.getfirst('pattern_id')

    # source_id becomes part of a file path, so allow only a plain name (no "../")
    if not source_id or secure_filename(source_id) != source_id or not pattern_id:
        print_error("400 Bad Request", "A valid source_id and pattern_id are required.")
        return

    file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format("cart." + source_id))

    try:
        df = pd.read_csv(file_path, sep="\t")
    except FileNotFoundError:
        print_error("404 Not Found", "This genecart is either not found, or is unweighted (every gene has the same weight)")
        return

    if pattern_id not in df.columns:
        print_error("404 Not Found", f"Pattern '{pattern_id}' was not found in this gene list.")
        return

    df.sort_values(by=pattern_id, ascending=False, inplace=True)

    # Col 0 is uniq ID, col 1 is gene symbol (but the column name may vary).
    gene_list = df.iloc[:,1].tolist()
    weight_list = df[pattern_id].tolist()

    result = []
    for idx, gene in enumerate(gene_list):
        weight = weight_list[idx]
        result.append({'gene':gene, "weight":weight})

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
