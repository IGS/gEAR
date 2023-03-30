#!/opt/bin/python3

"""
When a user uploads "three-tab" file  bundle it contains the following files:

  - expression.tab
  - genes.tab
  - observations.tab

It's expected that the values in the 1st column of the observations tab file
match positionally with the values in the row headers of the expression.tab
file.  If not, you just get this error from scanpy:

  Exception: Observation IDs from 'expressions' and 'observations' files are not the same.

Doesn't help you find where the issue is.  This script reports where they first disagree.
"""


import argparse
import itertools
import os

def main():
    parser = argparse.ArgumentParser( description="Checks that observation IDs match in the expression tab file")
    parser.add_argument('-o', '--observations_file', type=str, required=True, help='Path to an input observations tab file' )
    parser.add_argument('-e', '--expression_file', type=str, required=True, help='Path to an input expression file' )
    args = parser.parse_args()

    obs_ids = get_observation_ids(args.observations_file)
    
    for line in open(args.expression_file):
        line = line.rstrip()
        obs_columns = line.split("\t")[1:]
        # only the first line matters
        break

    print("Got {0} obs IDs from {1}".format(len(obs_ids), args.observations_file))
    print("Got {0} obs IDs from {1}".format(len(obs_columns), args.expression_file))

    if set(obs_ids) != set(obs_columns):
            raise Exception("Observation IDs from 'expressions' and 'observations' files are not the same.")
    
    col_num = 0
    
    for (of_id, ef_id) in zip(obs_ids, obs_columns):
        col_num += 1

        if of_id != ef_id:
            print("Obs column mismatch at row/column {0}. Expression file has ({1}) but observation file has ({2})".format(col_num, ef_id, of_id))
            sys.exit(1)

    print("All observation values seemed to match")

def get_observation_ids(infile):
    ids = list()
    line_count = 0

    for line in open(infile):
        line_count += 1

        if line_count == 1:
            continue
        
        ids.append(line.split('\t')[0])

    return ids

if __name__ == '__main__':
    main()
