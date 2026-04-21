## if you change this:
#   docker compose down
#   docker compose up -d
## Can check logs with:
#   docker logs -f gear-jupyterhu

import os
import sys
sys.path.insert(0, "/srv/jupyterhub")

from dockerspawner import DockerSpawner

HOST_JUPYTERHUB_ROOT = os.environ["HOST_JUPYTERHUB_ROOT"]
HOST_USERHOMES_ROOT = os.path.join(HOST_JUPYTERHUB_ROOT, "userhomes")
HOST_DATASETS_ROOT = os.environ["HOST_DATASETS_ROOT"]

c = get_config()

# -----------------------------------------------------------------------------
# Core Hub settings
# -----------------------------------------------------------------------------

c.JupyterHub.bind_url = "http://:8000"
c.JupyterHub.hub_ip = "0.0.0.0"

# Persist Hub db inside /data (mounted from ./state/hub-data)
c.JupyterHub.db_url = "sqlite:////data/jupyterhub.sqlite"

# Required for encrypted auth_state storage
c.Authenticator.enable_auth_state = True

# gEAR decides who is allowed and the launch token proves that
c.Authenticator.allow_all = True

# Use the custom launch-token authenticator
c.JupyterHub.authenticator_class = "gear_auth.GearLaunchTokenAuthenticator"

# -----------------------------------------------------------------------------
# DockerSpawner settings
# -----------------------------------------------------------------------------

c.JupyterHub.spawner_class = DockerSpawner
c.DockerSpawner.cmd = ["start-singleuser.py"]

# Remove containers when they stop
c.DockerSpawner.remove = False

# Use local Docker default bridge network unless you later define another one
c.DockerSpawner.network_name = "jupyterhub_default"
c.DockerSpawner.use_internal_ip = True

# Directory inside spawned notebook containers
c.DockerSpawner.notebook_dir = "/home/jovyan"

# Default image
c.DockerSpawner.image = "gear-notebook:r"

# Persistent home directories on host
# Host path is relative to where the Hub container sees it:
# ./userhomes is mounted into the Hub container at /srv/jupyterhub/userhomes
c.DockerSpawner.volumes = {}

# Run notebooks as a non-root user where possible
c.DockerSpawner.environment = {
    "CHOWN_HOME": "yes",
    "CHOWN_HOME_OPTS": "-R",
}

# Default resource limits
c.Spawner.cpu_limit = 2
c.Spawner.mem_limit = "16G"

# -----------------------------------------------------------------------------
# User-selectable profiles
# -----------------------------------------------------------------------------

c.Spawner.profile_list = [
    {
        "display_name": "Python (Scanpy)",
        "default": True,
        "spawner_override": {
            "image": "gear-notebook:py",
            "cpu_limit": 2,
            "mem_limit": "16G",
        },
    },
    {
        "display_name": "R (SingleCellExperiment / Seurat)",
        "spawner_override": {
            "image": "gear-notebook:r",
            "cpu_limit": 2,
            "mem_limit": "16G",
        },
    },
    {
        "display_name": "Big Memory (Python)",
        "spawner_override": {
            "image": "gear-notebook:py",
            "cpu_limit": 4,
            "mem_limit": "32G",
        },
    },
]

# -----------------------------------------------------------------------------
# Spawn hook: mount only authorized datasets
# -----------------------------------------------------------------------------

async def gear_pre_spawn_hook(spawner: DockerSpawner):
    auth_state = await spawner.user.get_auth_state()
    if not auth_state:
        raise RuntimeError("No auth_state available for user")

    datasets = auth_state.get("gear_datasets", [])
    selected_dataset = auth_state.get("gear_selected_dataset")
    notebook_env = auth_state.get("gear_notebook_env", "python")

    if notebook_env == "r":
        spawner.image = "gear-notebook:r"
    else:
        spawner.image = "gear-notebook:py"

    username = spawner.user.name
    user_home_host = os.path.join(HOST_USERHOMES_ROOT, username)

    os.makedirs(user_home_host, exist_ok=True)
    os.chown(user_home_host, 1000, 100)

    volumes = {
        user_home_host: {
            "bind": "/home/jovyan",
            "mode": "rw",
        }
    }

    if len(datasets) > 25:
        raise RuntimeError("Too many datasets requested for one session")

    mounted_filenames = []

    for host_path in datasets:
        allowed_prefix = os.path.join(HOST_DATASETS_ROOT, "")
        if not host_path.startswith(allowed_prefix):
            raise RuntimeError(f"Invalid dataset path: {host_path}")

        if not os.path.isfile(host_path):
            raise RuntimeError(f"Dataset file not found: {host_path}")

        basename = os.path.basename(host_path)
        container_path = f"/data/{basename}"

        volumes[host_path] = {
            "bind": container_path,
            "mode": "ro",
        }
        mounted_filenames.append(basename)

    spawner.volumes = volumes

    # Add convenience environment variables
    env = dict(spawner.environment or {})
    env["GEAR_DATASET_FILES"] = ":".join(mounted_filenames)

    if selected_dataset:
        env["GEAR_SELECTED_DATASET"] = os.path.basename(selected_dataset)

    env["GEAR_USERNAME"] = spawner.user.name
    spawner.environment = env

c.Spawner.pre_spawn_hook = gear_pre_spawn_hook

# -----------------------------------------------------------------------------
# Idle culler
# -----------------------------------------------------------------------------

c.JupyterHub.services = [
    {
        "name": "cull-idle",
        "admin": True,
        "command": [
            "python",
            "-m",
            "jupyterhub_idle_culler",
            "--timeout=3600",
            "--cull-every=300",
        ],
    }
]
