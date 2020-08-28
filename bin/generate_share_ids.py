#!/usr/bin/env python3

"""
Probably want to do this first:

SELECT id FROM dataset WHERE share_id = ''
INTO OUTFILE '/var/lib/mysql-files/generate_share_ids.data'
FIELDS TERMINATED BY '\t'
ENCLOSED BY ''
LINES TERMINATED BY '\n';


"""

import uuid

for line in open('generate_share_ids.data'):
    line = line.rstrip()
    share_id = str(uuid.uuid4())
    
    print("UPDATE dataset SET share_id = '{0}' WHERE id = '{1}';".format(share_id, line))
