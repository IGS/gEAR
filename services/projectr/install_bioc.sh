#!/bin/bash

set -x

Rmaj=R-4
Rver="${Rmaj}.5.2"

current_dir=$(pwd)

# Download R source file
curl -s -L -f -o /tmp/${Rver}.tar.gz  http://lib.stat.cmu.edu/R/CRAN/src/base/${Rmaj}/${Rver}.tar.gz || exit 1

# Verify the file was downloaded and is valid
tar -tzf /tmp/${Rver}.tar.gz > /dev/null || exit 1

# Extract
tar xzf /tmp/${Rver}.tar.gz -C /opt || exit 1

cd /opt/${Rver}
/opt/${Rver}/configure --with-readline=no --enable-R-shlib --enable-BLAS-shlib --with-x=no || exit 1
make -j$(nproc) || exit 1
make install || exit 1

# Run Rscript
cd ${current_dir}
Rscript --vanilla install_bioc.R || exit 1

# done