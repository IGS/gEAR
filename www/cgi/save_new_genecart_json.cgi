#!/opt/bin/python3

"""
Saves a GeneCart object (from genecart.js) serialized as
formData and submitted to this script.  This assumes we're creating
a NEW GeneCart object

This script is appropriate to call when passing JSON directly.
"""

import json
import sys
from pathlib import Path

TWO_LEVELS_UP = 2
abs_path_gear = Path(__file__).resolve().parents[TWO_LEVELS_UP]
abs_path_lib = abs_path_gear.joinpath('lib')
sys.path.insert(0, str(abs_path_lib))
import geardb

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

def main():

    # NOTE: Not going to add weighted gene and cart info, but may add subclasses in the future
    gc = geardb.GeneCart()
    data = json.load(sys.stdin)
    gc.update_from_json(data)

    # Weighted carts also get saved to file
    if gc.gctype == 'weighted-list':
        import scanpy as sc

        file_ext = ".tab"
        source_file_path = CARTS_BASE_DIR.joinpath("cart.{0}{1}".format(gc.share_id, file_ext))
        h5dest_file_path = CARTS_BASE_DIR.joinpath("cart.{0}.h5ad".format(gc.share_id))

        try:
            # Read from JSON properties instead of from the cart object
            # Write the tab file out
            with open(source_file_path, 'wb') as sfh:
                headers = ["dataRowNames"]
                headers.extend(data["weight_labels"])
                sfh.write(('\t'.join(headers) + '\n').encode())

                for gene in data["genes"]:
                    row = [gene["id"]]
                    weights = [str(w) for w in gene["weights"]]
                    row.extend(weights)
                    sfh.write(('\t'.join(row) + '\n').encode())
        except AttributeError as e:
            print("Status: 500 Internal Server Error")
            print("Content-Type: application/json;charset=utf-8\n")
            print("GeneCart object is not a weighted list")
            print(str(e))
            sys.exit()

        # Write the h5ad file out
        adata = sc.read_csv(source_file_path, delimiter="\t", first_column_names=True).transpose()
        adata.write(filename=h5dest_file_path)

    gc.save()

    result = { 'id': gc.id }
    print("Content-Type: application/json\n\n")
    print(json.dumps(result))


if __name__ == '__main__':
    main()
