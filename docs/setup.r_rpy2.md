# Setting up R, projectR, and the rpy2 module

The base version of R installed on Ubuntu Bionic (18.04) is not a high enough version to install projectR. So we need to install manually.

## Prerequesites to install via apt-get

`sudo apt-get install gfortran libbz2-dev libcurl4-openssl-dev liblzma-dev libpcre3 libpcre3-dev libgomp1 libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev`

## Installing R

Run this to install R, Bioconductor, and projectR:

```text
cd <gEAR\_git\_root>/services/projectr/
sudo sh ./install_bioc.sh
```

To ensure R's shared libraries are found create a file "libR.conf" in /etc/ld.so.conf.d and add the following contents:

```text
# libR default configuation
/usr/local/lib/R/lib
```

Then run `sudo ldconfig` to cache the shared libraries.  You can confirm shared libaries with `sudo ldconfig -v`

## Installing rpy2

`<python\_bin>/pip3 install rpy2==3.5.1`

Later versions seem to have an error that is not an issue with this version
