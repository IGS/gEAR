variable "TAG" {
  default = "latest"
}

variable "REGISTRY" {
  default = "adkinsrs"
}

function "date" {
    params = []
    result = "${formatdate("YYYY-MM-DD", timestamp())}" # later, convert to use formattimestamp
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
    "anndata_upload_consumer"
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

target "python-base" {
    context = "."
    dockerfile = "Dockerfile.python"
    tags = ["${REGISTRY}/gear-python-base:${TAG}", "${REGISTRY}/gear-python-base:${date()}"]
    platforms = ["linux/amd64", "linux/arm64"]
}

target "r-base" {
    context = "."
    dockerfile = "Dockerfile.r"
    tags = ["${REGISTRY}/gear-r-base:${TAG}", "${REGISTRY}/gear-r-base:${date()}"]
    platforms = ["linux/amd64", "linux/arm64"]
}

target "listener-python-base" {
    context = ".."
    dockerfile = "listeners/Dockerfile.python_base"
    tags = ["${REGISTRY}/gear-listener-python-base:${TAG}", "${REGISTRY}/gear-listener-python-base:${date()}"]
    platforms = ["linux/amd64", "linux/arm64"]
}

target "web" {
  context    = "."
  dockerfile = "Dockerfile"
  tags       = ["${REGISTRY}/umgear:${TAG}", "${REGISTRY}/umgear:${date()}"]
  platforms  = ["linux/amd64", "linux/arm64"]
}

target "panel" {
  context    = "../services/spatial"
  dockerfile = "Dockerfile"
  tags       = ["${REGISTRY}/spatial_panel_app:${TAG}", "${REGISTRY}/spatial_panel_app:${date()}"]
  platforms  = ["linux/amd64", "linux/arm64"]
}

target "gosling_upload_consumer" {
  context    = ".."
  dockerfile = "listeners/Dockerfile.gosling_upload"
  tags       = ["${REGISTRY}/gear_gosling_upload_consumer:${TAG}", "${REGISTRY}/gear_gosling_upload_consumer:${date()}"]
  platforms  = ["linux/amd64", "linux/arm64"]
}

target "anndata_upload_consumer" {
  context    = ".."
  dockerfile = "listeners/Dockerfile.anndata_upload"
  tags       = ["${REGISTRY}/gear_anndata_upload_consumer:${TAG}", "${REGISTRY}/gear_anndata_upload_consumer:${date()}"]
  platforms  = ["linux/amd64", "linux/arm64"]
}

# NOTE: Not currently implemented, so it is not in the groupings.
target "projectr_consumer" {
  context    = ".."
  dockerfile = "listeners/Dockerfile.projectr"
  tags       = ["${REGISTRY}/gear_projectr_consumer:${TAG}", "${REGISTRY}/gear_projectr_consumer:${date()}"]
  platforms  = ["linux/amd64", "linux/arm64"]
}