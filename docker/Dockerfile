############################################################
# Dockerfile to build gEAR web application
############################################################

# Stage 1: Build Python with --enable-shared
FROM ubuntu:jammy AS python-builder

ARG PYTHON_MAJOR_VERSION=3.14
ARG PYTHON_FULL_VERSION=${PYTHON_MAJOR_VERSION}.4
ARG PYTHON_URL=https://www.python.org/ftp/python/${PYTHON_FULL_VERSION}/Python-${PYTHON_FULL_VERSION}.tar.xz

# Install build dependencies
RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  apache2-dev \
  libffi-dev \
  libsqlite3-dev \
  libreadline-dev \
  libncursesw5-dev \
  libssl-dev \
  tk-dev \
  libgdbm-dev \
  libc6-dev \
  liblzma-dev \
  libbz2-dev \
  zlib1g-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libblosc-dev \
  libhdf5-dev \
  hdf5-helpers \
  hdf5-tools \
  libblas-dev \
  liblapack-dev \
  libcairo2-dev \
  libxml2-dev \
  libwebp-dev \
  ca-certificates \
  wget \
  gcc \
  g++ \
  make \
  pkg-config \
  git \
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# Build Python with shared library support
# --enable-shared is needed due to the Apache mod_wsgi module needing to link against the Python shared library.
WORKDIR /tmp
RUN wget $PYTHON_URL \
  && tar -xf Python-${PYTHON_FULL_VERSION}.tar.xz \
  && cd Python-${PYTHON_FULL_VERSION} \
  # LDFLAGS, CPPFLAGS, and PKG_CONFIG_PATH are needed to ensure the build process can find the lzma library and headers,
  # which are required for Python's lzma module.  Without these flags, the build process may fail to find lzma downstream.
  && ./configure --prefix=/opt/Python-${PYTHON_FULL_VERSION} --enable-optimizations --enable-shared \
    LDFLAGS="-L/usr/lib/x86_64-linux-gnu" \
    CPPFLAGS="-I/usr/include" \
    PKG_CONFIG_PATH="/usr/lib/x86_64-linux-gnu/pkgconfig" \
  && make install \
  && cd /tmp \
  && rm -rf Python-${PYTHON_FULL_VERSION}*

# Upgrade pip and install requirements from file
# This set of modules will be used by the main Python install and scripts that call it (like CGI ones)
# numpy is installed before requirements.txt because of dependency requirement issues
# NOTE: Installing custom version of diffxpy that is based on the latest commit on the main branch (at the time)
WORKDIR /opt/Python-${PYTHON_FULL_VERSION}/bin
COPY ./requirements.txt /opt/gEAR/docker/requirements.txt

ENV LD_LIBRARY_PATH="/opt/Python-${PYTHON_FULL_VERSION}/lib"

# Disable user-installs to guarantee packages go into /opt/Python...
ENV PIP_USER=0

# Install stuff
RUN ./python3 -m pip install --upgrade pip setuptools wheel \
  && ./python3 -m pip install -r /opt/gEAR/docker/requirements.txt --prefer-binary \
  && ./python3 -m pip uninstall dask-expr -y

# Cleanup
RUN chmod 777 /opt/Python-${PYTHON_FULL_VERSION}/lib/python${PYTHON_MAJOR_VERSION}/site-packages/scanpy/preprocessing/__pycache__ 2>/dev/null || true \
  && rm -rf ~/.cache/pip


# Stage 2: Install R packages using shell scirpts
FROM ubuntu:jammy AS r-builder

# originally used a bioconductor base image, but ran into an issue where
# the base image was using ubuntu:24.04 (instead of 22:04)
# which could cause gcc compiler issues in stage 3

# 1. Install required system dependencies for compiling R and packages
RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  curl \
  wget \
  build-essential \
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
  ca-certificates \
  tzdata \
  git \
  unzip \
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# 2. Install Bioconductor
WORKDIR /opt/gEAR/docker
COPY ./install_bioc.sh ./install_bioc.R ./
RUN chmod +x install_bioc.sh \
  && ./install_bioc.sh

# 3. Install other packages
COPY ./install_packages.R ./
RUN Rscript --vanilla install_packages.R || exit 1


# Stage 3: Final runtime image
FROM ubuntu:jammy AS final

ARG TARGETARCH
# Some package building issues when attempting to use 3.14
ARG PYTHON_MAJOR_VERSION=3.14
ARG PYTHON_FULL_VERSION=${PYTHON_MAJOR_VERSION}.4
ARG GEAR_URL=https://www.github.com/IGS/gEAR.git
ARG GEAR_BRANCH=devel

LABEL author="Shaun Adkins <sadkins@som.umaryland.edu>"
LABEL maintainer="Shaun Adkins <sadkins@som.umaryland.edu>"

# Install runtime packages
RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  ca-certificates \
  curl \
  make \
  gcc \
  g++ \
  rsync \
  wget \
  git \
  vim \
  unzip \
  llvm \
  apache2 \
  php \
  libffi-dev \
  libpcre2-dev \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libicu-dev \
  libdeflate-dev \
  libssl3 \
  libapache2-mod-php \
  php-gd \
  libffi-dev \
  libpcre2-dev \
  fonts-roboto \
  ttf-mscorefonts-installer \
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

ENV LLVM_CONFIG=/usr/bin/llvm-config-14

# Copy compiled Python from builder stage
COPY --from=python-builder /opt/Python-${PYTHON_FULL_VERSION} /opt/Python-${PYTHON_FULL_VERSION}

# Copy compiled R from r-builder stage
COPY --from=r-builder /usr/local/lib/R /usr/local/lib/R

# Link Python and shared library
RUN mkdir -p /opt/bin \
  && ln -s /opt/Python-${PYTHON_FULL_VERSION}/bin/python3 /opt/bin/python3 \
  && ln -s /opt/Python-${PYTHON_FULL_VERSION}/lib/libpython${PYTHON_MAJOR_VERSION}.so.1.0 /usr/lib/

# Tell the system linker about the new Python libs
RUN echo "/opt/Python-${PYTHON_FULL_VERSION}/lib" > /etc/ld.so.conf.d/python314.conf \
    && ldconfig

# LD_LIBRARY_PATH needs R library path for rpy2 compatibility and Python library path too
ENV LD_LIBRARY_PATH="/usr/local/lib/R/lib:/opt/Python-${PYTHON_FULL_VERSION}/lib"
ENV R_HOME="/usr/local/lib/R"
ENV PATH="/opt/Python-${PYTHON_FULL_VERSION}/bin:/usr/local/lib/R/bin:${PATH}"

# Force GCC 11 to use a newer C standard to understand R 4.5.2 headers
ENV CFLAGS="-std=gnu2x"

# Install rpy2 now that R has installed
RUN /opt/bin/python3 -m pip install rpy2==3.6.7

# Pre-compile all Python bytecode during the build so Apache doesn't try to write it at runtime
#RUN /opt/Python-${PYTHON_FULL_VERSION}/bin/python3 -m compileall -q /opt/Python-${PYTHON_FULL_VERSION}/lib/python${PYTHON_MAJOR_VERSION}/site-packages/
RUN find /opt/Python-${PYTHON_FULL_VERSION}/lib/python${PYTHON_MAJOR_VERSION}/site-packages/ -name __pycache__ -exec chmod 777 {} \;

# Get gEAR code.  Currently cannot do a git clone by normal means due to the repo being "private"
WORKDIR /opt
RUN git clone -b ${GEAR_BRANCH} --single-branch ${GEAR_URL} \
  && rm -rf /var/www \
  && ln -s /opt/gEAR/www/ /var/

COPY ./gear.ini.docker /opt/gEAR/gear.ini

# Apache things
WORKDIR /opt/gEAR
COPY ./apache2.conf /etc/apache2/apache2.conf
COPY ./wsgi.conf /etc/apache2/mods-available/wsgi.conf

# The wsgi.load .so file name is different for different architectures,
# so we have separate ones for each architecture.
COPY ./wsgi.${TARGETARCH}.load /etc/apache2/mods-available/wsgi.load
COPY ./000-default.conf /etc/apache2/sites-available/000-default.conf
COPY ./umgear.conf /etc/apache2/sites-available/umgear.conf
COPY ./php.ini /etc/php/8.1/apache2/php.ini
RUN a2dissite 000-default.conf
RUN a2ensite umgear.conf
RUN a2enmod rewrite cgi wsgi proxy proxy_http proxy_wstunnel include headers
EXPOSE 80

CMD ["/usr/sbin/apachectl", "-D", "FOREGROUND"]
