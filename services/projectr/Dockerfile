
# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM python:3.10

# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED=True

RUN apt-get -qq update \
    && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
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
    # Cleanup
    && apt -qq clean autoclean \
    && apt -qq autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# ProjectR install.  Must be done before rpy2 is installed
# NOTE: Speed this up... so slow for Linux install.  Maybe build on actual linux hardware, or do a multi-stage docker build
COPY ./install_bioc.sh /tmp/
COPY ./install_bioc.R /tmp/
WORKDIR /tmp
RUN /tmp/install_bioc.sh

# Install production dependencies.
COPY ./requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

ENV LD_LIBRARY_PATH="/usr/local/lib/R/lib:$LD_LIBRARY_PATH"

# Copy local code to the container image.
ENV APP_HOME=/app
WORKDIR $APP_HOME
COPY . ./

# Run the web service on container startup. Here we use the gunicorn
# webserver, with one worker process and 8 threads.
# For environments with multiple CPU cores, increase the number of workers
# to be equal to the cores available.
# Timeout is set to 0 to disable the timeouts of the workers to allow Cloud Run to handle instance scaling.
CMD exec gunicorn --bind :$PORT --workers 4 --threads 8 --timeout 0 main:app