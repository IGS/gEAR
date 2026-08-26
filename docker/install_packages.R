#!/usr/bin/env Rscript --vanilla

options(warn = 2) # Treat warnings as errors (i.e. missing Ubuntu modules)

# Tell pak to prefer binary packages
# NOTE: Unset pkg.sysreqs if on server, as apt update will fail otherwise
options(pkg.sysreqs = TRUE)
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))
#options(Ncpus = 1)

Sys.setenv(PKG_INCLUDE_LINKINGTO = "TRUE")  # Set this environment variable to avoid issues with packages that link to other packages

install.packages("pak", ask=FALSE, repos="https://cloud.r-project.org/")
library(pak)

# Install in smaller batches to avoid memory/subprocess issues
packages_list <- list(c(
        'reticulate'
        , 'ctlab/fgsea'
        , 'genesofeve/projectR@d3dd79e2b14172a9561059d58462c97f0a78d4c8'
        , 'stuart-lab/signac'
        , 'httpuv'
        , 'hdf5r'
        , 'Seurat'
        , 'bioc::rhdf5'
        , 'bioc::anndataR'
        , 'bioc::biomaRt'
        , "CHuanSite/SJD"
        )
)

for (pkg_batch in packages_list) {
    message("\n=== Installing batch: ", paste(pkg_batch, collapse=", "), " ===")
    tryCatch( {
        pak::pkg_install(pkg_batch)
    }, error = function(e) {
        message("Error: ", conditionMessage(e))
        if (requireNamespace("rlang", quietly = TRUE)) {
            print(rlang::last_trace(drop = FALSE))
        }
        quit(status = 1, save = "no")
    })
}