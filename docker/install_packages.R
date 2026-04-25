#!/usr/bin/env Rscript --vanilla

options(warn = 2) # Treat warnings as errors (i.e. missing Ubuntu modules)

# Bioconductor should already be installed.
# Load libraries
library(remotes)    # for install_version

tryCatch( {
    remotes::install_version("reticulate", version="1.46.0", repos="https://cloud.r-project.org/", ask=FALSE, dependencies=NA) # Sanity check with rpy2
    remotes::install_github("ctlab/fgsea")   # needed for projectR
    remotes::install_github("genesofeve/projectR@d3dd79e2b14172a9561059d58462c97f0a78d4c8")  # version 1.23.2
    BiocManager::install("biomaRt", ask=FALSE) # version 2.60.0
    remotes::install_github("CHuanSite/SJD")
    }, error = function(e) {
        message("Error: ", e$message)
        quit(status = 1, save = "no")
    }
)
