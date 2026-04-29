# Setting up R, projectR, and the rpy2 module

Rather than using the OS' package manager to install R, we'll use the steps/script below to do it so we can more carefully control versioning.

## Prerequesites to install via apt

```bash
sudo apt-get -qq update
sudo DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  gfortran \
  libreadline-dev \
  xorg-dev \
  libbz2-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libfontconfig1-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libwebp-dev \
  libgit2-dev \
  libuv1-dev \
  tzdata \
sudo apt -qq clean autoclean
sudo apt -qq autoremove -y
sudo rm -rf /var/lib/apt/lists/*
```

## Installing R

Run this to install R, Bioconductor, and projectR:

```bash
cd <gEAR\_git\_root>/services/projectr/
sudo sh ./install_bioc.sh
sudo Rscript --vanilla install_packages.R || exit 1

export R_HOME="/usr/local/lib/R"
export LD_LIBRARY_PATH="/usr/local/lib/R/lib:${LD_LIBRARY_PATH}"

```

To ensure R's shared libraries are found create a file "libR.conf" in /etc/ld.so.conf.d and add the following contents:

```bash
# libR default configuation
/usr/local/lib/R/lib
```

Then run `sudo ldconfig` to cache the shared libraries.  You can confirm shared libaries with `sudo ldconfig -v`

## Installing rpy2

This will be done in the Python installation doc steps.
