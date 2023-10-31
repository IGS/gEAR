#!/bin/bash

Rmaj=R-4
Rver="${Rmaj}.3.1"

current_dir=$(pwd)

# Install and build R (Using 'apt-get install' on Ubuntu Trusty installs version 3.0.2 of R)
curl http://lib.stat.cmu.edu/R/CRAN/src/base/${Rmaj}/${Rver}.tar.gz | tar -C /opt -zx
cd /opt/${Rver}
/opt/${Rver}/configure --with-readline=no --enable-R-shlib --enable-BLAS-shlib --with-x=no || exit 1
make || exit 1
make install || exit 1

# Run Rscript
cd ${current_dir}
Rscript --vanilla install_bioc.R || exit 1

# done