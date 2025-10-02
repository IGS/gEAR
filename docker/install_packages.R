#!/usr/bin/env Rscript --vanilla

# Bioconductor should already be installed.
# Load libraries
library(remotes)    # for install_version

install_version("devtools", version="2.4.5", repos="https://cloud.r-project.org/", ask=FALSE)
library(devtools)   # for install_github

tryCatch( {
    install.packages("reticulate", repos = "https://cloud.r-project.org/", ask = FALSE) # Sanity check with rpy2
    install_github("ctlab/fgsea")   # needed for projectR
    install_github("genesofeve/projectR@d3dd79e2b14172a9561059d58462c97f0a78d4c8")  # version 1.23.2
    BiocManager::install("biomaRt", ask=FALSE) # version 2.60.0
    install_github("CHuanSite/SJD")
    }, error = function(e) {
        message("Error: ", e$message)
        quit(status = 1, save = "no")
    }
)
