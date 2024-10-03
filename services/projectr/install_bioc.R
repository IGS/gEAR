#!/usr/bin/env Rscript --vanilla

install.packages(c("BiocManager", "devtools"), dependencies=TRUE, repos="http://lib.stat.cmu.edu/R/CRAN/")
BiocManager::install(version = "3.19")  # required for R 4.4.0
BiocManager::install(c("genesofeve/projectR", "biomaRt"), ask=FALSE)
library(devtools); install_github("CHuanSite/SJD")