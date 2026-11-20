variable "TAG" {
  default = "latest"
}

variable "REGISTRY" {
  default = "adkinsrs"
}

# Grab date from the passed in envvar set in the terminal
# DATE=$(date +%Y-%m-%d) docker buildx bake --allow=fs.read=.. <group>
variable "DATE" {
    default = "devel"
}

group "all" {
  targets = [
    "python-base",
    "r-base",
    "listener-python-base",
    "web",
    "panel",
    "gosling_upload_consumer",
    "anndata_upload_consumer"
  ]
}

# These are the final images, where only the codebase has changed.
group "default" {
  targets = [
    "web",
    "panel",
    "gosling_upload_consumer",
    "anndata_upload_consumer",
    "spatial_upload_consumer"
  ]
}

# These are intermediate images that are used in the build process of the final images.
group "intermediate" {
    targets = [
        "python-base",
        "r-base",
        "listener-python-base",
    ]
}

# Define a reusable base configuration target
target "_common" {
    platforms = ["linux/amd64", "linux/arm64"]
    output    = ["type=registry"]  # auto-pushes to dockerhub, but MUST be pulled locally.
}

target "python-base" {
    inherits = ["_common"]
    context = "."
    dockerfile = "Dockerfile.python"
    tags = ["${REGISTRY}/gear-python-base:${TAG}", "${REGISTRY}/gear-python-base:${DATE}"]
}

target "r-base" {
    inherits = ["_common"]
    context = "."
    dockerfile = "Dockerfile.r"
    tags = ["${REGISTRY}/gear-r-base:${TAG}", "${REGISTRY}/gear-r-base:${DATE}"]
}

target "listener-python-base" {
    inherits = ["_common"]
    context = ".."
    dockerfile = "listeners/Dockerfile.python_base"
    tags = ["${REGISTRY}/gear-listener-python-base:${TAG}", "${REGISTRY}/gear-listener-python-base:${DATE}"]
}

target "web" {
    inherits = ["_common"]
    context    = "."
    dockerfile = "Dockerfile"
    tags       = ["${REGISTRY}/umgear:${TAG}", "${REGISTRY}/umgear:${DATE}"]
}

target "panel" {
    inherits = ["_common"]
    context    = "../services/spatial"
    dockerfile = "Dockerfile"
    tags       = ["${REGISTRY}/spatial_panel_app:${TAG}", "${REGISTRY}/spatial_panel_app:${DATE}"]
}

target "gosling_upload_consumer" {
    inherits = ["_common"]
    context    = ".."
    dockerfile = "listeners/Dockerfile.gosling_upload"
    tags       = ["${REGISTRY}/gear_gosling_upload_consumer:${TAG}", "${REGISTRY}/gear_gosling_upload_consumer:${DATE}"]
}

target "anndata_upload_consumer" {
    inherits = ["_common"]
    context    = ".."
    dockerfile = "listeners/Dockerfile.anndata_upload"
    tags       = ["${REGISTRY}/gear_anndata_upload_consumer:${TAG}", "${REGISTRY}/gear_anndata_upload_consumer:${DATE}"]
}

target "spatial_upload_consumer" {
    inherits = ["_common"]
    context    = ".."
    dockerfile = "listeners/Dockerfile.spatial_upload"
    tags       = ["${REGISTRY}/gear_spatial_upload_consumer:${TAG}", "${REGISTRY}/gear_spatial_upload_consumer:${DATE}"]
}

# NOTE: Not currently implemented, so it is not in the groupings.
target "projectr_consumer" {
    inherits = ["_common"]
    context    = ".."
    dockerfile = "listeners/Dockerfile.projectr"
    tags       = ["${REGISTRY}/gear_projectr_consumer:${TAG}", "${REGISTRY}/gear_projectr_consumer:${DATE}"]
}

target "projectr_cloud_run" {
    platforms = ["linux/amd64"] # Only push unix image, as Google Cloud Run can only use linux/amd64 images.
    output    = ["type=registry"]
    context = "../services/projectr"
    dockerfile = "Dockerfile"
    tags = ["us-east1-docker.pkg.dev/gear-154704/cloud-run-source-deploy/projectr_service:${TAG}",
    "us-east1-docker.pkg.dev/nemo-analytics/cloud-run-source-deploy/projectr_service:${TAG}",
    "us-east1-docker.pkg.dev/inflammation-gear/cloud-run-source-deploy/projectr_service:${TAG}"
    ]

}