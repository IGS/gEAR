############################################################
# Dockerfile to build gEAR web application
############################################################

# Stage 1: Build Python with --enable-shared
FROM ubuntu:jammy AS python-builder

ARG PYTHON_VERSION=3.14.4
ARG PYTHON_URL=https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz

# Install build dependencies
RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  ca-certificates \
  libffi-dev \
  libsqlite3-dev \
  libreadline-dev \
  libncursesw5-dev \
  libssl-dev \
  tk-dev \
  libgdbm-dev \
  libc6-dev \
  libbz2-dev \
  zlib1g-dev \
  wget \
  build-essential \
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# Build Python with shared library support
# --enable-shared is needed due to the Apache mod_wsgi module needing to link against the Python shared library.
WORKDIR /tmp
RUN wget $PYTHON_URL \
  && tar -xf Python-${PYTHON_VERSION}.tar.xz \
  && cd Python-${PYTHON_VERSION} \
  && ./configure --prefix=/opt/Python-${PYTHON_VERSION} --enable-optimizations --enable-shared \
  && make install \
  && cd /tmp \
  && rm -rf Python-${PYTHON_VERSION}*

# Stage 2: Install R packages (using Bioconductor base image)
FROM bioconductor/bioconductor_docker:3.21-r-4.5.2 AS r-builder

# Install additional runtime dependencies if needed
RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  ca-certificates \
  curl \
  # to build rpy2
  python3-pip \
  build-essential \
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# Copy install scripts
COPY ./install_packages.R /tmp/

# Install R packages
RUN /usr/local/bin/Rscript --vanilla /tmp/install_packages.R || exit 1

# Pre-build rpy2 wheel
RUN pip3 wheel --no-cache-dir --no-deps --wheel-dir /tmp/wheels rpy2==3.6.7

# Stage 3: Final runtime image
FROM ubuntu:jammy

ARG TARGETARCH
# Some package building issues when attempting to use 3.14
ARG PYTHON_MAJOR_VERSION=3.14
ARG PYTHON_FULL_VERSION=${PYTHON_MAJOR_VERSION}.4
ARG CMAKE_VERSION=3.19.8
ARG CMAKE_URL=https://github.com/Kitware/CMake/archive/refs/tags/v${CMAKE_VERSION}.tar.gz
ARG GEAR_URL=https://www.github.com/IGS/gEAR.git
ARG GEAR_BRANCH=devel

# For pybigwig dependency
ENV PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1

LABEL maintainer="Shaun Adkins <sadkins@som.umaryland.edu>"

# Install runtime packages
RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  gcc \
  g++ \
  make \
  libffi-dev \
  libsqlite3-dev \
  libhdf5-dev \
  libreadline-dev \
  libncursesw5-dev \
  libssl-dev \
  tk-dev \
  libgdbm-dev \
  libc6-dev \
  libbz2-dev \
  hdf5-helpers \
  hdf5-tools \
  zlib1g-dev \
  libblas-dev \
  liblapack-dev \
  libcairo2-dev \
  libxml2-dev \
  libhdf5-serial-dev \
  libwebp-dev \
  ca-certificates \
  curl \
  rsync \
  wget \
  git \
  vim \
  unzip \
  llvm \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libgit2-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libblosc-dev \
  apache2-dev \
  apache2 \
  php \
  libapache2-mod-php \
  php-gd \
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

ENV LLVM_CONFIG=/usr/bin/llvm-config-14

# Copy compiled Python from builder stage
COPY --from=python-builder /opt/Python-${PYTHON_FULL_VERSION} /opt/Python-${PYTHON_FULL_VERSION}

# Copy compiled R from r-builder stage
COPY --from=r-builder /usr/local/lib/R /usr/local/lib/R

# Copy pre-built rpy2 wheel
COPY --from=r-builder /tmp/wheels /tmp/wheels


# Link Python and shared library
RUN mkdir -p /opt/bin \
  && ln -s /opt/Python-${PYTHON_FULL_VERSION}/bin/python3 /opt/bin/python3 \
  && ln -s /opt/Python-${PYTHON_FULL_VERSION}/lib/libpython3.14.so.1.0 /usr/lib/

# LD_LIBRARY_PATH needs R library path for rpy2 compatibility
ENV LD_LIBRARY_PATH="/usr/local/lib/R/lib"
ENV R_HOME="/usr/local/lib/R"

# Get gEAR code.  Currently cannot do a git clone by normal means due to the repo being "private"
WORKDIR /opt
RUN git clone -b ${GEAR_BRANCH} --single-branch ${GEAR_URL} \
  && rm -rf /var/www \
  && ln -s /opt/gEAR/www/ /var/

COPY ./gear.ini.docker /opt/gEAR/gear.ini

# Upgrade pip and install requirements from file
# This set of modules will be used by the main Python install and scripts that call it (like CGI ones)
# numpy is installed before requirements.txt because of dependency requirement issues
# NOTE: Installing custom version of diffxpy that is based on the latest commit on the main branch (at the time)
WORKDIR /opt/Python-${PYTHON_FULL_VERSION}/bin
COPY ./requirements.txt /opt/gEAR/docker/requirements.txt
RUN ./pip3 install --upgrade pip setuptools wheel \
  # Ran into "error: resolution-too-deep" so used `pip freeze` from a production env to create requirements.full.txt
  # Source: https://medium.com/@ebenyemiriam17/error-resolution-too-deep-how-i-fixed-my-pip-dependency-resolution-error-9cb3471344f9
  && ./pip3 install -r /opt/gEAR/docker/requirements.txt  --prefer-binary \
  # Uninstall dask-expr to avoid issues with spatialdata (https://github.com/scverse/spatialdata/pull/570)
  && ./pip3 uninstall dask-expr -y \
  && (chmod 777 /opt/Python-${PYTHON_FULL_VERSION}/lib/python${PYTHON_MAJOR_VERSION}/site-packages/scanpy/preprocessing/__pycache__ || true) \
  && rm -rf ~/.cache/pip \
  && find /opt/Python-${PYTHON_FULL_VERSION} -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

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