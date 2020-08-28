#!/usr/bin/env python3

import io
import re
import sys

replacements = dict()
replacements_str = """\
anatomy_id_3 Cochlear_vascular_endothelium
anatomy_id_4 Cochlear_mesenchyme
anatomy_id_5 Cochlear_neurons
anatomy_id_6 Cochlear_epithelial_cells
anatomy_id_7 Cochlear_epithelial_nonsensory
anatomy_id_8 Cochlear_epithelial_hair_cells
anatomy_id_9 Vestibular_vascular_endothelium
anatomy_id_10 Vestibular_mesenchyme
anatomy_id_11 Vestibular_neurons
anatomy_id_12 Vestibular_epithelial_cells
anatomy_id_13 Vestibular_epithelial_nonsensory
anatomy_id_14 Vestibular_epithelial_hair_cells
"""

rfh = io.StringIO(replacements_str)
for d in rfh:
    anatomy, name = d.strip().split(' ')
    replacements[anatomy] = name

for line in open(sys.argv[1]):
    for k in replacements:
        line = re.sub("\"{0}\"".format(k),"\"{0}\"".format(replacements[k]), line.rstrip())
    print(line)
