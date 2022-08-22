#!/opt/bin/python3

"""

"""

import cgi
import json
from pathlib import Path

abs_path_www = Path(__file__).resolve().parents[1] # web-root dir
CARTS_BASE_DIR = abs_path_www.joinpath("carts")

def main():
    form = cgi.FieldStorage()
    source_id = form.getvalue('source_id')  # Root of the file name (minus extension)

    result = []

    print('Content-Type: application/json\n\n')

    # Handle unweighted genecarts which are not saved to tabfile.
    if not "cart." in source_id:
        result.append({"label":"unweighted"})
        print(json.dumps(result))
        return

    # TODO: Consider loading from h5ad instead of tab if it exists
    file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format(source_id))

    for line in open(file_path):
        line = line.rstrip()
        # on the first line, the file should be all pattern names after the first column
        cols = line.split("\t")

        # Col 0 is uniq ID, col 1 is gene symbol.
        for col in cols[2:]:
            result.append({'label': col})

        # we only care about the first line
        break

    print(json.dumps(result))

if __name__ == '__main__':
    main()
