#!/opt/bin/python3

"""
Uses the EpiViz gene lookup utility written by Jayaram to provide information about 
gene's structural annotation in mouse.  This prevent us from having to sync our
annotation version with that maintained by EpiViz.

Gene symbols are searched case-insensitively.

Output example:

{
  requestId: 0,
  type: "response",
  data: [
     {
        probe: "MGI:5477246",
        gene: "Gm26752",
        seqName: "10",
        start: 3364073,
        end: 3367135
     }
  ]
}

These are the sort of entries returned by the Epiviz utility.

  {'end': 72622829, 'gene': 'Rfx7', 'seqName': '9', 'start': 72532240, 'probe': 'MGI:2442675'}
  {'end': 72470756, 'gene': 'Rfx7', 'seqName': 'chr9', 'start': 72380046}

Here we ignore any entries with a 'probe' key, and filter for only the first
exact match 'gene' entry.
"""

import cgi
import json
import requests
import sys

def main():
    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    gene_symbol = form.getvalue('gene_symbol')
    
    #r = requests.get("http://epiviz-dev.cbcb.umd.edu/data2/main.php?action=search&q={0}".format(gene_symbol))
    r = requests.get("https://epiviz.umgear.org/app/hrp/igs/data/main.php?action=search&q={0}".format(gene_symbol))

    # By default the EpiViz gene lookup utility returns all matches where the
    #  query is a substring.  Filter these to only get the right exact match
    jdata = r.json()

    for row in jdata['data']:
        if 'probe' not in row:
            if gene_symbol.lower() == row['gene'].lower():
                print(json.dumps(row))

if __name__ == '__main__':
    main()
