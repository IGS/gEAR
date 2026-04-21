#!/bin/bash

export JUPYTERHUB_CRYPT_KEY="$(openssl rand -hex 32)"
export GEAR_LAUNCH_SECRET="$(openssl rand -hex 32)"
export HUB_PUBLIC_BASE_URL="http://localhost"
