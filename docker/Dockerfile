############################################################
# Dockerfile to build gEAR web application
############################################################

FROM ubuntu:bionic

LABEL maintainer="Shaun Adkins <sadkins@som.umaryland.edu>"

ARG PYTHON_VERSION=3.7.3
ARG PYTHON_URL=https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz

# To assign $git_user run `docker build --build-arg git_user=<user> --build-arg git_pass=<pass> .`
# Avoiding using my personal username and password for this
ARG GEAR_BRANCH=main
ARG git_user
ARG git_pass
ARG GEAR_URL=https://${git_user}:${git_pass}@github.com/IGS/gEAR.git

# Install packages to Ubuntu
RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt-get -qq install -y --no-install-recommends \
  libffi-dev \
  libsqlite-dev \
  libsqlite3-dev \
  libhdf5-dev \
  libreadline-gplv2-dev \
  libncursesw5-dev \
  libssl-dev python3-tk \
  tk-dev \
  libgdbm-dev \
  libc6-dev \
  libbz2-dev \
  hdf5-helpers \
  hdf5-tools \
  libhdf5-dev \
  zlib1g-dev \
  libblas-dev \
  liblapack-dev \
  libxml2-dev \
  cmake \
  curl \
  wget \
  git \
  vim \
  unzip \
  software-properties-common \
  build-essential \
  python-dev \
  openssh-server \
  # Apache2 and PHP
  apache2-dev \
  apache2 \
  php \
  libapache2-mod-php \
  # Cleanup
  && apt-get -qq clean autoclean \
  && apt-get -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# Install Python (with configuration options needed for gEAR)
# TODO: test using a apt install instead, now that virtualenv is not needed
WORKDIR /tmp
RUN wget $PYTHON_URL \
  && tar -xf Python-${PYTHON_VERSION}.tar.xz \
  && cd Python-${PYTHON_VERSION} \
  && ./configure --prefix=/opt/Python-${PYTHON_VERSION} --enable-optimizations --enable-shared \
#  && make \
  && make install \
  && cd /tmp \
  && rm -rf  Python-${PYTHON_VERSION}* \
  && mkdir -p /opt/bin \
  && ln -s /opt/Python-${PYTHON_VERSION}/bin/python3 /opt/bin/python3

# Get gEAR code.  Currently cannot do a git clone by normal means due to the repo being "private"
WORKDIR /opt
RUN git clone -b ${GEAR_BRANCH} --single-branch ${GEAR_URL} \
  && rm -rf /var/www \
  && ln -s /opt/gEAR/www/ /var/

COPY ./gear.ini.docker /opt/gEAR/gear.ini

# Upgrade pip and install requirements from file
# This set of modules will be used by the main Python install and scripts that call it (like CGI ones)
# NOTE: Installing custom version of diffxpy that I forked and fixed a function to work
WORKDIR /opt/Python-${PYTHON_VERSION}/bin
COPY ./requirements.txt /opt/gEAR/docker/requirements.txt
RUN ln -s /opt/Python-${PYTHON_VERSION}/lib/libpython3.7m.so.1.0 /usr/lib/ \
  && ./pip3 install --upgrade pip \
  && ./pip3 install -r /opt/gEAR/docker/requirements.txt \
  && /opt/Python-${PYTHON_VERSION}/bin/python3 -m pip install git+https://github.com/adkinsrs/diffxpy.git@b2ebeb0fb7c6c215d51264cd258edf9d013ff021 \
  && apt-get -qq update \
  && apt-get -qq install -y --no-install-recommends libapache2-mod-wsgi \
  && (chmod 777 /opt/Python-${PYTHON_VERSION}/lib/python${PYTHON_VERSION}/site-packages/scanpy/preprocessing/__pycache__ || true)

# Apache things
WORKDIR /opt/gEAR
COPY ./apache2.conf /etc/apache2/apache2.conf
COPY ./wsgi.conf /etc/apache2/mods-available/wsgi.conf
COPY ./wsgi.load /etc/apache2/mods-available/wsgi.load
COPY ./000-default.conf /etc/apache2/sites-available/000-default.conf
COPY ./umgear.conf /etc/apache2/sites-available/umgear.conf
COPY ./php.ini /etc/php/7.2/apache2/php.ini
RUN a2dissite 000-default.conf
RUN a2ensite umgear.conf
RUN a2enmod rewrite cgi wsgi proxy
EXPOSE 80

CMD /usr/sbin/apachectl -D FOREGROUND
