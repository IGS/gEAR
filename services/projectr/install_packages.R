#!/usr/bin/env Rscript --vanilla

# Bioconductor should already be installed.

# Load libraries
library(devtools)
library(remotes)

tryCatch( {
    install.packages("reticulate", repos = "https://cloud.r-project.org/", ask = FALSE) # Sanity check with rpy2
    install_github("ctlab/fgsea")   # needed for projectR
    install_github("genesofeve/projectR")
    BiocManager::install("biomaRt", ask=FALSE) # version 2.60.0
    install_github("CHuanSite/SJD")
    }, error = function(e) {
        message("Error: ", e$message)
        quit(status = 1, save = "no")
    }
)
