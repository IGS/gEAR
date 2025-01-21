# Use the official lightweight Python image.
# https://hub.docker.com/_/python
FROM python:3.12

# Allow statements and log messages to immediately appear in the Knative logs
ENV PYTHONUNBUFFERED=True

RUN apt-get -qq update \
  && DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  gcc \
  libblosc-dev \
  # Cleanup
  && apt -qq clean autoclean \
  && apt -qq autoremove -y \
  && rm -rf /var/lib/apt/lists/*

# Copy local code to the container image.
ENV APP_HOME=/app
WORKDIR $APP_HOME
COPY . ./

# Install production dependencies.
RUN pip install -r requirements.txt

# Run the web service on container startup.
CMD ["panel", "serve", "panel_app.py", "--address", "0.0.0.0", "--port", "8081", "--num-threads", "4", "--allow-websocket-origin", "*"]

# Note: I set this to 8081 so I can use as a standalone image, where 8080 is used by the main app