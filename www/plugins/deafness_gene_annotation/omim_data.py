#!/usr/bin/env python3

import json

data = dict()

for line in open('./omim_data.tab'):
    line = line.rstrip()
    cols = line.split("\t")

    genes = cols[2].split(',')

    for gene in genes:
        gene = gene.replace(' ', '').lower()

        if gene not in data:
            data[gene] = {
                'phenotypes': [],
                'links': [{'label': 'OMIM', 'url': "https://omim.org/entry/{0}".format(cols[4])}]
            }

        data[gene]['phenotypes'].append(cols[8])

print(json.dumps(data, indent=3))
