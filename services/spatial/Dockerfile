# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM python:3.12

# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED=True

# NOTE: Installing vim so I can debug the container
RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  gcc \
  libblosc-dev \
  vim \
  # Cleanup
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# Install production dependencies.
COPY ./requirements.txt ./
RUN pip install -r requirements.txt
# Uninstall dask-expr to avoid issues with spatialdata (https://github.com/scverse/spatialdata/pull/570)
RUN pip uninstall dask-expr -y

# Copy local code to the container image.
# App home is reflecting the relative levels within the gEAR repo
ENV APP_HOME=/gEAR/services/spatial
WORKDIR $APP_HOME
COPY . ./

EXPOSE 5006

# Run the web service on container startup.
# NOTE: num-threads is set in the app, and will throw an error if also set here
CMD ["panel", "serve", "panel_app.py", "panel_app_expanded.py", "--address", "0.0.0.0", "--port", "5006", "--allow-websocket-origin", "*", "--global-loading-spinner"]

# For the future.  Is added to systemd service file
#"--num-procs", "4",

# Things I would like to test -
#, "--reuse-sessions", "--rest-session-info"]