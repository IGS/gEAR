#!/bin/bash

## Exports all the SQL needed to replicate a dataset onto another system, useful 
#  for development.  You still need to migrate files like the H5AD file on your
#  own.

ACCOUNT_INFO_FILE=$HOME/.my.cnf

# they needed to pass a dataset ID
DATASET_ID=$1

if [ "$DATASET_ID" == "" ]; then
    echo "ERROR: Pass an dataset ID"
    exit 1
fi

# make sure they have a credentials file
#  https://stackoverflow.com/a/9293090
if [ ! -f "$ACCOUNT_INFO_FILE" ]; then
    echo "ERROR: Expected credentials file $ACCOUNT_INFO_FILE to exist"
    exit 1
fi

## do them in order so they can be copy/pasted
mysqldump --compact -t -h localhost gear_portal dataset --where="id='$DATASET_ID'"
mysqldump --compact -t -h localhost gear_portal dataset_display --where="dataset_id='$DATASET_ID'"

#mysqldump --compact -t -h localhost gear_portal dataset_preference --where="dataset_id='$DATASET_ID'"
#mysqldump --compact -t -h localhost gear_portal dataset_shares --where="dataset_id='$DATASET_ID'"
