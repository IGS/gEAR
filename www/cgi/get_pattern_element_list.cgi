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
    file_path = Path(CARTS_BASE_DIR).joinpath("{}.tab".format(source_id))

    for line in open(file_path):
        line = line.rstrip()
        # on the first line, the file should be all pattern names after the first column
        cols = line.split("\t")

        for col in cols[1:]:
            result.append({'label': col})

        # we only care about the first line
        break

    print('Content-Type: application/json\n\n')
    print(json.dumps(result))

if __name__ == '__main__':
    main()
