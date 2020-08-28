#!/usr/bin/env python3

"""
Prints a summary of the use of a dataset based on apache logs.  

Expectations of log entries:

Analyser loads
70.169.55.242 - - [23/Jan/2020:14:45:15 +0000] "GET /analyze_dataset.html?dataset_id=c4f16a12-9e98-47be-4335-b8321282919e 
HTTP/1.1"

Comparison tool
134.192.135.254 - - [23/Jan/2020:18:31:51 +0000] "GET /compare_datasets.html?dataset_id=c4f16a12-9e98-47be-4335-b8321282919e HTTP/1.1" 200 8097

Front-page views
72.213.153.240 - - [23/Jan/2020:19:13:15 +0000] "GET /api/plot/c4f16a12-9e98-47be-4335-b8321282919e

"""

import os
import re
import sys
import gzip

dataset_id = sys.argv[1]
APACHE_LOG_BASE = '/var/log/apache2'

if not dataset_id:
    print("\nUsage: report_dataset_usage_stats.py [dataset_id]")
    sys.exit(1)

view_users = dict()
view_count = 0
    
analyzer_users = dict()
analyzer_count = 0

comparison_users = dict()
comparison_count = 0

for file in os.listdir(APACHE_LOG_BASE):
    m = re.match(".+access.log.*", file)
    if m:
        filepath = "{0}/{1}".format(APACHE_LOG_BASE, file)

        if filepath.endswith('.gz'):
            fh = gzip.open(filepath, 'rb')
            is_compressed = True
        else:
            fh = open(filepath, 'r')
            is_compressed = False

        for line in fh:
            if is_compressed:
                line = line.decode()

            # front page views
            m = re.search("^(\S+).+GET \/api\/plot\/(\S+)\/", line)
            if m:
                ip_address = m.group(1)
                this_dataset_id = m.group(2)

                if dataset_id == this_dataset_id:
                    view_count += 1
                    view_users[ip_address] = 1

            # match analyzer load
            m = re.search("^(\S+).+GET \/analyze_dataset.html.dataset_id=(\S+)", line)
            if m:
                ip_address = m.group(1)
                this_dataset_id = m.group(2)

                if dataset_id == this_dataset_id:
                    analyzer_count += 1
                    analyzer_users[ip_address] = 1

            # match comparison tool loads
            m = re.search("^(\S+).+GET \/compare_datasets.html.dataset_id=(\S+)", line)
            if m:
                ip_address = m.group(1)
                this_dataset_id = m.group(2)

                if dataset_id == this_dataset_id:
                    comparison_count += 1
                    comparison_users[ip_address] = 1

print("Viewed on front page by {0} times by {1} users".format(view_count, len(view_users)))
print("Single-cell workbench used {0} times by {1} users".format(analyzer_count, len(analyzer_users)))
print("Comparison used {0} times by {1} users".format(comparison_count, len(comparison_users)))
