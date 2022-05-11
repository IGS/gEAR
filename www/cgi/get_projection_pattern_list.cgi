#!/opt/bin/python3

"""
Used by projection.html, this script gets a list of the available analysis pattern files.
Returns as a list
"""

import cgi, json
from pathlib import Path

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

def main():
    #form = cgi.FieldStorage()

    result = []

    for p_file in CARTS_BASE_DIR.iterdir():
        name = p_file.name    # get basname
        if not name.endswith(".tab"):   # skip non-tab files, like .gitignore
            continue
        # NOTE: Weighted gene carts may be retrieved here in the future as well.
        if name.startswith("cart."):    # skip weighted gene carts, as those are handled separately.
            continue
        source_id = name.replace('.tab', '')
        title = source_id.replace('ROWmeta_DIMRED_', '')
        result.append({"id":source_id, "title":title})

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
