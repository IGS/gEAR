#!/opt/bin/python3

"""
Lists the available DominoSignal output files in www/cache/domino_signal and, if a
filename is passed, the selectable parameter values within that file and the dataset details
from its matching .json file (e.g. domino_network_p21.reformatted.json), if there is one.

Data structure returned:

{
    success: 1,
    files: ['domino_network_p21.reformatted.parquet', ...],
    options: {                       # only when filename is passed
        tfs: [...],
        ligands: [...],
        receptors: [...],
        sending_clusters: [...],
        receiving_clusters: [...]
    },
    dataset_info: {                  # only when filename is passed; null if there is no matching .json file
        title, author, year, doi, pmid, gear_url, description   # each a string, or null if blank
    }
}

"""

import cgi
import json
import os
import sys

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDOUT, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

from gear import dominosignal as ds


def main():
    form = cgi.FieldStorage()
    filename = form.getvalue('filename')

    result = {'success': 1, 'files': ds.list_domino_files()}

    if filename:
        result['options'] = ds.get_domino_options(ds.load_domino_df(filename))
        result['dataset_info'] = ds.load_dataset_info(filename)

    return result


if __name__ == '__main__':
    try:
        result = main()
    except ds.DominoFileError as e:
        result = {'success': -1, 'message': str(e)}
    except Exception as e:
        result = {'success': -1, 'message': "Error reading DominoSignal file: {0}".format(e)}

    sys.stdout = original_stdout
    print('Content-Type: application/json\n\n')
    print(json.dumps(result))
