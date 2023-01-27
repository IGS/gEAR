#!/opt/bin/python3

"""
For a given weighted genecart, return all genes and their weights (sorted descending) for the specified pattern
"""

import cgi
import json
import pandas as pd
from pathlib import Path

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

def main():
    form = cgi.FieldStorage()
    source_id = form.getvalue('source_id')  # Root of the file name (minus extension)
    pattern_id = form.getvalue('pattern_id')

    print('Content-Type: application/json\n\n')

    try:
        file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format("cart." + source_id))

        result = []

        df = pd.read_csv(file_path, sep="\t")
        df.sort_values(by=pattern_id, ascending=False, inplace=True)

        # Col 0 is uniq ID, col 1 is gene symbol (but the column name may vary).
        gene_list = df.iloc[:,1].tolist()
        weight_list = df[pattern_id].tolist()

        for idx, gene in enumerate(gene_list):
            weight = weight_list[idx]
            result.append({'gene':gene, "weight":weight})

        print(json.dumps(result))
    except FileNotFoundError:
        print("This genecart is either not found, or is unweighted (every gene has the same weight)")

if __name__ == '__main__':
    main()
