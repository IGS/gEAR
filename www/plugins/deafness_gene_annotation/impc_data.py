#!/usr/bin/env python3

import json

genes = dict()

for line in open('./impc_data.tab'):
    line = line.rstrip()
    cols = line.split("\t")

    genes[cols[0].lower()] = {
        'on_hover': 'IMPC - ' + cols[1],
        'links': [{'label': 'IMPC', 'url': cols[2]}]
    }

print(json.dumps(genes, indent=3))
