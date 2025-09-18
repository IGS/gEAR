#!/usr/bin/env Rscript --vanilla

# Install required packages
tryCatch( {
    install.packages(c("BiocManager", "remotes"), dependencies=TRUE, repos="http://lib.stat.cmu.edu/R/CRAN/")
    BiocManager::install(version = "3.19", ask=FALSE)  # required for R 4.4.0
    }, error = function(e) {
        message("Error: ", e$message)
        quit(status = 1, save = "no")
    }
)
# Individual packages are installed in a separate script so that the lengthy installation process can be cached
