############################################################
# Dockerfile to build gEAR web application
############################################################

FROM ubuntu:jammy
# ubuntu 22.04

LABEL maintainer="Shaun Adkins <sadkins@som.umaryland.edu>"

ARG CMAKE_VERSION=3.19.8
ARG CMAKE_URL=https://github.com/Kitware/CMake/archive/refs/tags/v${CMAKE_VERSION}.tar.gz

ARG PYTHON_VERSION=3.10.4
ARG PYTHON_URL=https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz

ARG GEAR_URL=https://www.github.com/IGS/gEAR.git
ARG GEAR_BRANCH=devel

# Install packages to Ubuntu
RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  libffi-dev \
  libsqlite-dev \
  libsqlite3-dev \
  libhdf5-dev \
  libreadline-dev \
  libncursesw5-dev \
  libssl-dev \
  tk-dev \
  libgdbm-dev \
  libc6-dev \
  libbz2-dev \
  lzma \
  liblzma-dev \
  hdf5-helpers \
  hdf5-tools \
  libhdf5-dev \
  zlib1g-dev \
  libblas-dev \
  liblapack-dev \
  libcairo2-dev \
  libxml2-dev \
  libhdf5-serial-dev \
  libwebp-dev \
  cmake \
  curl \
  rsync \
  wget \
  git \
  vim \
  unzip \
  software-properties-common \
  build-essential \
  gcc \
  python3-dev \
  python3-tk \
  openssh-server \
  llvm \
  # Required for R
	gfortran \
  # Required for rpy2
  r-cran-rjava \
  # Required for R-package devtools (which is required for SJD)
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libgit2-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  # Required for Panel (though we use Panel in a separate image)
  libblosc-dev \
  # Apache2 and PHP
  libapache2-mod-wsgi-py3 \
  apache2-dev \
  apache2 \
  php \
  libapache2-mod-php \
  php-gd \
  # Cleanup
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

ENV LLVM_CONFIG=/usr/bin/llvm-config-9

# Install CMake
# Ubuntu Bionic package is 3.10 which is too outdated for louvain-igraph install
# The Kitware PPA repo for apt does not go earlier than 3.20, and anything later than 3.20 also breaks
# So need to build from source instead

WORKDIR /tmp
RUN wget $CMAKE_URL \
  && tar -xzf v${CMAKE_VERSION}.tar.gz \
  && cd CMake-${CMAKE_VERSION} \
  && ./bootstrap \
  && make install \
  && rm -rf CMake-${CMAKE_VERSION} v${CMAKE_VERSION}.tar.gz

# ProjectR install.  Must be done before rpy2 is installed
# NOTE: Speed this up... so slow for Linux install.  Maybe build on actual linux hardware, or do a multi-stage docker build
COPY ./install_bioc.sh /tmp/
COPY ./install_bioc.R /tmp/
COPY ./install_packages.R /tmp/
WORKDIR /tmp
RUN /tmp/install_bioc.sh
RUN Rscript --vanilla install_packages.R || exit 1

# Install Python (with configuration options needed for gEAR)
# TODO: test using a apt install instead, now that virtualenv is not needed
WORKDIR /tmp
RUN wget $PYTHON_URL \
  && tar -xf Python-${PYTHON_VERSION}.tar.xz \
  && cd Python-${PYTHON_VERSION} \
  && ./configure --prefix=/opt/Python-${PYTHON_VERSION} --enable-optimizations --enable-shared \
  && make install \
  && cd /tmp \
  && rm -rf  Python-${PYTHON_VERSION}* \
  && mkdir -p /opt/bin \
  && ln -s /opt/Python-${PYTHON_VERSION}/bin/python3 /opt/bin/python3 \
  && ln -s /opt/Python-${PYTHON_VERSION}/lib/libpython3.10.so.1.0 /usr/lib/

# Get gEAR code.  Currently cannot do a git clone by normal means due to the repo being "private"
WORKDIR /opt
RUN git clone -b ${GEAR_BRANCH} --single-branch ${GEAR_URL} \
  && rm -rf /var/www \
  && ln -s /opt/gEAR/www/ /var/

COPY ./gear.ini.docker /opt/gEAR/gear.ini

# Add some executables
WORKDIR /opt/gEAR
RUN mkdir -p src \
  && rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/hubClone /opt/gEAR/src/ \
  && rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/hubCheck /opt/gEAR/src/ \
  && rsync -aP hgdownload.soe.ucsc.edu::genome/admin/exe/linux.x86_64/bigBedToBed /opt/gEAR/src/


# Upgrade pip and install requirements from file
# This set of modules will be used by the main Python install and scripts that call it (like CGI ones)
# numpy is installed before requirements.txt because of dependency requirement issues
# NOTE: Installing custom version of diffxpy that is based on the latest commit on the main branch (at the time)
WORKDIR /opt/Python-${PYTHON_VERSION}/bin
COPY ./requirements.full.txt /opt/gEAR/docker/requirements.full.txt
RUN ./pip3 install --upgrade pip \
  # wheel is needed for the installation of biopython
  && ./pip3 install wheel \
  # Ran into "error: resolution-too-deep" so used `pip freeze` from a production env to create requirements.full.txt
  # Source: https://medium.com/@ebenyemiriam17/error-resolution-too-deep-how-i-fixed-my-pip-dependency-resolution-error-9cb3471344f9
  && ./pip3 install -r /opt/gEAR/docker/requirements.full.txt \
  # Uninstall dask-expr to avoid issues with spatialdata (https://github.com/scverse/spatialdata/pull/570)
  && ./pip3 uninstall dask-expr -y \
  && apt-get -qq update \
  && (chmod 777 /opt/Python-${PYTHON_VERSION}/lib/python${PYTHON_VERSION}/site-packages/scanpy/preprocessing/__pycache__ || true)

# Scikit-learn lib -> https://stackoverflow.com/a/70016416
#ENV LD_PRELOAD="/usr/lib/aarch64-linux-gnu/libgomp.so.1:/opt/Python-${PYTHON_VERSION}/lib/python3.10/site-packages/scikit_learn.libs/libgomp-d22c30c5.so.1.0.0"
ENV LD_LIBRARY_PATH="/usr/local/lib/R/lib:$LD_LIBRARY_PATH"

# Apache things
WORKDIR /opt/gEAR
COPY ./apache2.conf /etc/apache2/apache2.conf
COPY ./wsgi.conf /etc/apache2/mods-available/wsgi.conf
COPY ./wsgi.load /etc/apache2/mods-available/wsgi.load
COPY ./000-default.conf /etc/apache2/sites-available/000-default.conf
COPY ./umgear.conf /etc/apache2/sites-available/umgear.conf
COPY ./php.ini /etc/php/8.1/apache2/php.ini
RUN a2dissite 000-default.conf
RUN a2ensite umgear.conf
RUN a2enmod rewrite cgi wsgi proxy proxy_http proxy_wstunnel include headers
EXPOSE 80

CMD /usr/sbin/apachectl -D FOREGROUND