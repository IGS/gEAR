#!/bin/bash

Rmaj=R-4
Rver="${Rmaj}.4.0"

current_dir=$(pwd)

curl -s -L http://lib.stat.cmu.edu/R/CRAN/src/base/${Rmaj}/${Rver}.tar.gz | tar xzv -C /opt
cd /opt/${Rver}
/opt/${Rver}/configure --with-readline=no --enable-R-shlib --enable-BLAS-shlib --with-x=no || exit 1
make || exit 1
make install || exit 1

# Run Rscript
cd ${current_dir}
Rscript --vanilla install_bioc.R || exit 1

# done