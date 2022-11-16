#!/bin/bash

Rmaj=R-4
Rver="${Rmaj}.2.0"  # Consider moving up to v4

# Install and build R (Using 'apt-get install' on Ubuntu Trusty installs version 3.0.2 of R)
curl https://cran.r-project.org/src/base/${Rmaj}/${Rver}.tar.gz | tar -C /opt -zx
cd /opt/${Rver}
/opt/${Rver}/configure --with-readline=no --enable-R-shlib --enable-BLAS-shlib --with-x=no || exit 1
make || exit 1
make install || exit 1

apt-get -qq install -y --no-install-recommends r-base-dev

# install projectR via Bioconductor
echo "install.packages(c(\"BiocManager\"), repos=\"http://lib.stat.cmu.edu/R/CRAN/\")" | R --save --restore || exit 1
echo "BiocManager::install(c(\"projectR\"), ask=FALSE)" | R --save --restore || exit 1

# done