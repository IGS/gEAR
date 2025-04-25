#!/bin/bash

# This script is used to serve the panel apps in the background
# It is intended to be run as a systemd service with a load balancer

# Path to the panel executable
PANEL_PATH=$1
if [ -z "$PANEL_PATH" ]; then
    echo "Usage: $0 <path_to_panel_executable>"
    exit 1
fi

# Serve the panel apps
${PANEL_PATH} serve panel_app.py --port 5006 --address=0.0.0.0  --num-procs=4 --allow-websocket-origin="*" --global-loading-spinner &
${PANEL_PATH} serve panel_app.py --port 5007 --address=0.0.0.0  --num-procs=4 --allow-websocket-origin="*" --global-loading-spinner &
/${PANEL_PATH} serve panel_app.py --port 5008 --address=0.0.0.0  --num-procs=4 --allow-websocket-origin="*" --global-loading-spinner &

# Never stop
wait -n