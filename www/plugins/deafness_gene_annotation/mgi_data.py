#!/usr/bin/env python3

import json

genes = dict()

for line in open('./mgi_data.tab'):
    line = line.rstrip()
    cols = line.split("\t")

    genes[cols[0].lower()] = {
        'phenotypes': list(),
        'links': [{'label': 'MGI',
                   'url': "http://www.informatics.jax.org/marker/summary?nomen={0}&cm=&coordinate=&coordUnit=bp&startMarker=&endMarker=&go=&goVocab=goFunctionTerm&goVocab=goProcessTerm&goVocab=goComponentTerm&interpro=&phenotype=".format(cols[0])}]
    }

    if len(cols) == 3:
        genes[cols[0].lower()]['phenotypes'] = cols[2].split(' | ')

print(json.dumps(genes, indent=3))
