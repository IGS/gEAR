#!/usr/bin/env Rscript --vanilla

# Install required packages
install.packages(c("BiocManager", "devtools", "remotes"), dependencies=TRUE, repos="http://lib.stat.cmu.edu/R/CRAN/")
BiocManager::install(version = "3.19", ask=FALSE)  # required for R 4.4.0

# Individual packages are installed in a separate script so that the lengthy installation process can be cached
