#!/opt/bin/python3

"""
For a given weighted or unweighted genecart, return the pattern labels.

If genecart is weighted, return the top 10 up- and down-regulated genes for each pattern
"""

import cgi, sys
import json
import pandas as pd
from pathlib import Path

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

def main():
    form = cgi.FieldStorage()
    source_id = form.getvalue('source_id')  # Root of the file name (minus extension)
    scope = form.getvalue("scope", None)
    result = []

    print('Content-Type: application/json\n\n')

    # Handle unweighted genecarts which are not saved to tabfile.
    # NOTE: Does not check db for existence of unweighted genecart
    if scope == "unweighted-list":
        result.append({"label":"unweighted", "top_up":"", "top_down":"", "binary": True})
        print(json.dumps(result))
        return

    # TODO: Consider loading from h5ad instead of tab if it exists
    file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format("cart." + source_id))

    try:
        df = pd.read_csv(file_path, sep="\t")
    except FileNotFoundError as e:
        print(e, file=sys.stderr)
        print(json.dumps(result))
        return

    # Col 0 is uniq ID, col 1 is gene symbol (but the column name may vary).
    for col in df.columns[2:]:
        up_genes = df.nlargest(n=5, columns=[col]).iloc[:, 1].tolist()
        down_genes = df.nsmallest(n=5, columns=[col]).iloc[:, 1].tolist()

        # if there are no negative values, the down_genes should be empty
        # Example is with binary weights.
        if all([v >= 0 for v in df[col]]):
            down_genes = []

        # If all values are either 0 or 1, the pattern can be a binary weight
        binary = all([v in [0, 1] for v in df[col]])

        result.append({'label': col, "top_up":",".join(up_genes), "top_down":",".join(down_genes), "binary": binary})

    print(json.dumps(result))

if __name__ == '__main__':
    main()
