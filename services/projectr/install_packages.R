#!/usr/bin/env Rscript --vanilla

# Bioconductor should already be installed.
# Load libraries
library(remotes)    # for install_version

# Install devtools dependencies explicitly (these will have R 4.5-compatible binaries)
tryCatch({
    install.packages(c("usethis", "fs", "miniUI", "pak", "pkgdown", "pkgload", "profvis", "roxygen2", "testthat", "urlchecker"),
                     repos="https://cloud.r-project.org/", ask=FALSE)
}, error = function(e) {
    message("Error installing devtools dependencies: ", e$message)
    quit(status = 1, save = "no")
})

tryCatch({
    remotes::install_version("devtools", version="2.5.1", repos="https://cloud.r-project.org/", ask=FALSE, dependencies=FALSE)
    #install.packages("devtools", repos = "https://cloud.r-project.org/", ask = FALSE, dependencies = TRUE)
}, error = function(e) {
    message("Error installing devtools: ", e$message)
    quit(status = 1, save = "no")
})
library(devtools)   # for install_github

tryCatch( {
    remotes::install_version("reticulate", version="1.46.0", repos="https://cloud.r-project.org/", ask=FALSE, dependencies=TRUE) # Sanity check with rpy2
    #install.packages("reticulate", repos = "https://cloud.r-project.org/", ask = FALSE, dependencies = TRUE)
    install_github("ctlab/fgsea")   # needed for projectR
    install_github("genesofeve/projectR@d3dd79e2b14172a9561059d58462c97f0a78d4c8")  # version 1.23.2
    BiocManager::install("biomaRt", ask=FALSE) # version 2.60.0
    install_github("CHuanSite/SJD")
    }, error = function(e) {
        message("Error: ", e$message)
        quit(status = 1, save = "no")
    }
)
