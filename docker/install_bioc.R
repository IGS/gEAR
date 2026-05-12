#!/usr/bin/env Rscript --vanilla

# Install required packages
tryCatch( {
    install.packages(c("BiocManager", "remotes"), dependencies=NA, repos="http://lib.stat.cmu.edu/R/CRAN/")
    BiocManager::install(version = "3.21", ask=FALSE)
    }, error = function(e) {
        message("Error: ", e$message)
        quit(status = 1, save = "no")
    }
)
# Individual packages are installed in a separate script so that the lengthy installation process can be cached
