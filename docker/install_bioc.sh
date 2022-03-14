#!/bin/bash
# R version for ubuntu bionic is 3.4.4, so use old method of installation
echo "source(\"http://bioconductor.org/biocLite.R\"); biocLite(c(\"projectR\"), ask=FALSE)" | R --save --restore || exit 1

# done