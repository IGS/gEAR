# check=skip=JSONArgsRecommended

FROM python:3.14-slim

# Immediately flush logs
ENV PYTHONUNBUFFERED=True

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libgomp1 \
    libcurl4 \
    libxml2 \
    libdeflate0 \
    libtirpc3 \
    liblapack3 \
    libblas3 \
    libpng16-16 \
    && rm -rf /var/lib/apt/lists/*

# Copy pre-built artifacts from r-builder stage
COPY --from=adkinsrs/gear-r-base:latest /usr/lib/x86_64-linux-gnu/libicuuc.so.70 /usr/lib/x86_64-linux-gnu/
COPY --from=adkinsrs/gear-r-base:latest /usr/lib/x86_64-linux-gnu/libicui18n.so.70 /usr/lib/x86_64-linux-gnu/
COPY --from=adkinsrs/gear-r-base:latest /usr/lib/x86_64-linux-gnu/libicudata.so.70 /usr/lib/x86_64-linux-gnu/
COPY --from=adkinsrs/gear-r-base:latest /usr/local/lib/R /usr/local/lib/R
COPY --from=adkinsrs/gear-r-base:latest /usr/local/bin/R /usr/local/bin/R
COPY --from=adkinsrs/gear-r-base:latest /usr/local/bin/Rscript /usr/local/bin/Rscript

# LD_LIBRARY_PATH needs R library path for rpy2 compatibility
RUN ldconfig
ENV PATH="/usr/local/bin:${PATH}"
ENV R_HOME="/usr/local/lib/R"
ENV LD_LIBRARY_PATH="/usr/local/lib/R/lib"

# Force GCC 11 to use a newer C standard to understand R 4.x headers
ENV CFLAGS="-std=gnu2x"

# Install production dependencies.
COPY ./requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt \
  && rm -rf ~/.cache/pip
RUN find /opt/Python-${PYTHON_FULL_VERSION} -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

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
CMD exec python -m gunicorn --bind :${PORT:-8080} --workers 3 --threads 3 --timeout 0 main:app