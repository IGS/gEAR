#!/usr/bin/env Rscript --vanilla

install.packages(c("BiocManager", "devtools"), dependencies=TRUE, repos="http://lib.stat.cmu.edu/R/CRAN/")
BiocManager::install(c("projectR", "biomaRt"), ask=FALSE)
library(devtools); install_github("CHuanSite/SJD")