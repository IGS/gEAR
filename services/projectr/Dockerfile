
# Stage 1: Install R packages (using Bioconductor base image)
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


# Stage 2: Final image
FROM python:3.14-slim

# Copy compiled R from r-builder stage
COPY --from=r-builder /usr/local/lib/R /usr/local/lib/R

# Copy pre-built rpy2 wheel
COPY --from=r-builder /tmp/wheels /tmp/wheels

# LD_LIBRARY_PATH needs R library path for rpy2 compatibility
ENV LD_LIBRARY_PATH="/usr/local/lib/R/lib"
ENV R_HOME="/usr/local/lib/R"

# Install production dependencies.
COPY ./requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt \
  && rm -rf ~/.cache/pip \
  && find /opt/Python-${PYTHON_FULL_VERSION} -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

# Copy local code to the container image.
ENV APP_HOME=/app
WORKDIR $APP_HOME
COPY . ./

# Run the web service on container startup. Here we use the gunicorn
# webserver, with one worker process and 8 threads.
# For environments with multiple CPU cores, increase the number of workers
# to be equal to the cores available.
# Timeout is set to 0 to disable the timeouts of the workers to allow Cloud Run to handle instance scaling.

# SAdkins - "2 * CPU + 1" formula to set Gunicorn workers and threads together
# https://medium.com/building-the-system/gunicorn-3-means-of-concurrency-efbb547674b7
CMD ["exec", "gunicorn", "--bind", ":%{PORT:-8080}", "--workers", "3", "--threads", "3", "--timeout", "0", "main:app"]