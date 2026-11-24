# syntax=docker/dockerfile:1
############################################################
# Dockerfile to build gEAR web application
############################################################

# Stage 3: Final runtime image
FROM ubuntu:jammy AS final

ARG TARGETARCH
# Some package building issues when attempting to use 3.14
ARG PYTHON_MAJOR_VERSION=3.14
ARG PYTHON_FULL_VERSION=${PYTHON_MAJOR_VERSION}.4
ARG GEAR_URL=https://github.com/IGS/gEAR.git
ARG GEAR_BRANCH=devel

LABEL author="Shaun Adkins <sadkins@som.umaryland.edu>"
LABEL maintainer="Shaun Adkins <sadkins@som.umaryland.edu>"

# Install runtime packages
RUN apt -qq update \
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
  mysql-client \
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
  libgfortran5 \
  libuv1 \
  libcairo2 \
  libhdf5-dev \
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# Cache the Roboto and Arial fonts
RUN fc-cache -f -v

ENV LLVM_CONFIG=/usr/bin/llvm-config-14

# Copy compiled Python from builder stage
COPY --from=adkinsrs/gear-python-base:latest /opt/Python-${PYTHON_FULL_VERSION} /opt/Python-${PYTHON_FULL_VERSION}

# Copy compiled R from r-builder stage
COPY --from=adkinsrs/gear-r-base:latest /usr/local/lib/R /usr/local/lib/R
COPY --from=adkinsrs/gear-r-base:latest /usr/local/bin/R /usr/local/bin/R
COPY --from=adkinsrs/gear-r-base:latest /usr/local/bin/Rscript /usr/local/bin/Rscript

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

# Force GCC 11 to use a newer C standard to understand R 4.x headers
ENV CFLAGS="-std=gnu2x"

# Install rpy2 now that R has installed
RUN /opt/bin/python3 -m pip install rpy2==3.6.7

# Pre-compile all Python bytecode during the build so Apache doesn't try to write it at runtime
#RUN /opt/Python-${PYTHON_FULL_VERSION}/bin/python3 -m compileall -q /opt/Python-${PYTHON_FULL_VERSION}/lib/python${PYTHON_MAJOR_VERSION}/site-packages/
RUN find /opt/Python-${PYTHON_FULL_VERSION}/lib/python${PYTHON_MAJOR_VERSION}/site-packages/ -name __pycache__ -exec chmod 777 {} \;

# Chromium for Plotly static image export (fig.to_image / write_image, e.g. the "Download image" menu item).
# kaleido 1.x does not include a browser. Its own installer (plotly_get_chrome / kaleido_get_chrome) is not used because:
#  - it installs into the installing user's home, which the Apache user (www-data) cannot see
#  - Google publishes no Linux arm64 build of the Chrome it downloads
# Playwright publishes Chromium for both amd64 and arm64, and --with-deps installs the system libraries it needs.
# The browser goes in a shared directory, linked to /usr/bin/chromium, one of the paths kaleido checks.
# Playwright itself is only needed to download the browser, so it goes in a throwaway venv.
# PLAYWRIGHT_BROWSERS_PATH is exported inside the RUN so the step also works when pasted into a shell
# (otherwise Playwright downloads to the current user's ~/.cache/ms-playwright).
ARG PLAYWRIGHT_VERSION=1.55.0
ARG CHROMIUM_DIR=/opt/ms-playwright
RUN export PLAYWRIGHT_BROWSERS_PATH=${CHROMIUM_DIR} \
  && /opt/bin/python3 -m venv /tmp/playwright-venv \
  && /tmp/playwright-venv/bin/pip install -q playwright==${PLAYWRIGHT_VERSION} \
  && apt -qq update \
  && /tmp/playwright-venv/bin/playwright install --with-deps --no-shell chromium \
  && CHROME_BIN=$(find ${CHROMIUM_DIR} -type f -name chrome -path '*/chromium-*' 2>/dev/null | head -n 1) \
  && { [ -n "${CHROME_BIN}" ] || { echo "ERROR: Chromium executable not found under ${CHROMIUM_DIR}" >&2; exit 1; }; } \
  && ln -sf "${CHROME_BIN}" /usr/bin/chromium \
  && chmod -R a+rX ${CHROMIUM_DIR} \
  && rm -rf /tmp/playwright-venv /var/lib/apt/lists/*

# Get gEAR code. ADD with a git source checks the branch's latest commit on every build and
# rebuilds this layer when the branch has moved. A "RUN git clone" layer would be reused from
# the build cache and keep serving old code.
# --keep-git-dir keeps .git so git commands still work inside the container.
ADD --keep-git-dir=true ${GEAR_URL}#${GEAR_BRANCH} /opt/gEAR
RUN rm -rf /var/www \
  && ln -s /opt/gEAR/www/ /var/

# gear.ini comes from the Docker template in the cloned repo, so no secrets are baked into the
# published image. Its [database] section already matches the "db" service in docker-compose.yml;
# the other placeholders (e-mail, GitHub token, etc.) stay unfilled.
# docker-compose.yml bind-mounts the host checkout over /opt/gEAR, which hides this copy,
# so the host checkout needs its own gear.ini (see docs/developer/setup/docker.md).
RUN cp /opt/gEAR/docker/gear.ini.docker.template /opt/gEAR/gear.ini

# Site preferences come from the template in the cloned repo, pointed at the local Docker site.
# domain_label "gear" matches the "portals" lists of the bundled study pages (www/landing/*/study.json).
# docker-compose.yml bind-mounts the host checkout over /opt/gEAR, which hides this copy,
# so the host checkout needs its own www/site_domain_prefs.json (see docs/developer/setup/docker.md).
# www/stats_archive/intro_stats.json is not made here; it needs the database (bin/generate_intro_stats.py)
# and the pages that read it hide those sections until it exists.
RUN /opt/bin/python3 -c "\
import json; \
prefs_path = '/opt/gEAR/www/site_domain_prefs.json'; \
prefs = json.load(open(prefs_path + '.template')); \
prefs.update(domain_label='gear', domain_short_display_label='gEAR', domain_url='http://localhost:8080'); \
json.dump(prefs, open(prefs_path, 'w'), indent=4)"

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
