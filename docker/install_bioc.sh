#!/bin/bash

echo "source(\"http://bioconductor.org/biocLite.R\"); biocLite(c(\"projectR\"), ask=FALSE)" | R --save --restore || exit 1

# done