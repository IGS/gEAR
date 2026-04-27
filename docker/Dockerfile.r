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
