############################################################
# Dockerfile to build gEAR web application
############################################################

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
  fontconfig \
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# Cache the Roboto and Arial fonts
RUN fc-cache -f -v

ENV LLVM_CONFIG=/usr/bin/llvm-config-14

# Copy compiled Python from builder stage
COPY --from=adkinsrs/gear-python-base:2026-04-27 /opt/Python-${PYTHON_FULL_VERSION} /opt/Python-${PYTHON_FULL_VERSION}

# Copy compiled R from r-builder stage
COPY --from=adkinsrs/gear-r-base:2026-04-27 /usr/local/lib/R /usr/local/lib/R
COPY --from=adkinsrs/gear-r-base:2026-04-27 /usr/local/bin/R /usr/local/bin/R
COPY --from=adkinsrs/gear-r-base:2026-04-27 /usr/local/bin/Rscript /usr/local/bin/Rscript

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
