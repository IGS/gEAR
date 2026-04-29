# Notes about the Docker setup

## Before building

* [Orbstack](https://orbstack.dev/) is recommended over Docker Desktop, but certainly optional. It is faster, lighter, and still uses the same Docker command-line tools.
  * Unfortunately Orbstack does not have Windows support.
* From the gEAR root, `cd docker`. All commands assume you are in this directory
* `cp docker-compose.yml.template docker-compose.yml`
  * Alternatively ask @adkinsrs for a docker-compose.yml file as it will be filled in. Otherwise fill in any values wrapped in brackets
  * @adkinsrs's file is hard-coded to his paths so be sure to change those.

## Acquiring the gEAR image

There are two options here.  The first method is significantly quicker, but there is a chance it may not be updated based on the latest Dockerfile instructions in Github (if @adkinsrs forgets to push the latest image up).

### Method 1: Pull image

Recently, docker images have switched to a multi-platform build where Docker will internally determine which image to pull based on your current platform. This requires you to ensure that the containerd store option is enabled on your copy of Docker Desktop.

Info on how to enable the containerd store

* https://docs.docker.com/engine/storage/containerd/

How to pull the image

* `docker pull adkinsrs/umgear:latest`.

Grab the gear.ini file:

* `docker run -it adkinsrs/umgear:latest`
* `docker ps` note the container ID
* `docker cp <container_id>:/opt/gEAR/gear.ini <gear_root>/gear.ini` which will copy the gear.ini file within the docker container outside.
* Feel free to stop the container, or leave it up so the docker-compose stack can use it.

This step is necessary if you are mounting your host gEAR code as a volume to the docker compose "web" service.  When you do this, the directory inside the container is replaced with your volume's code, and the gear.ini file within will be missing.

### Method 2: Build image

* Ensure you are in the "devel" branch of gEAR before building (`git checkout devel`)
* `cp gear.ini.docker.template gear.ini.docker`
  * Alternatively ask @adkinsrs for a gear.ini.docker file as it will be filled in. Otherwise fill in any values wrapped in brackets
* To build run `docker buildx build -t umgear:latest .`
  * Ensure the image name here (`umgear:latest`) is reflected in the docker-compose.yml file instead of `adkinsrs/umgear:latest`

### Method 2a: Build image with updated R or Python stuff

* This will use premade Python and R base Dockerfile images to save on build time. If you want to update the R or Python install, you need to do the following:
  * Update requirements.txt (for Python) or install_bioc.R, install_bioc.sh, or install_packages.R as needed for R
  * For R, run `docker buildx build -t gear-r-base -f Dockerfile.r .` (enjoy the bioconductor install slowness)
  * For Python, run  `docker buildx build -t gear-python-base -f Dockerfile.python .`
  * In the Dockerfile, change the `COPY --from` commands to point to your reviews r-base or python-base images
* To build gEAR, run `docker buildx build -t umgear:latest .`
  * Ensure the image name here (`umgear:latest`) is reflected in the docker-compose.yml file instead of `adkinsrs/umgear:latest`

In the build, the "gear.ini.docker" file will end up copied to "gear.ini" in the "/opt/gEAR" directory for the docker instance. However, if are using docker-compose and the gEAR directory is mounted into the "web" service, this can be overriden to a gear.ini from outside.  If you do not have a "gear.ini" file (only gear.ini.template), then ask @adkinsrs for one.

### The three Dockerfiles

Information about the three Dockerfiles found in `<gear_root>/docker`

#### Dockerfile.python (The Python Base)

This file is dedicated entirely to compiling Python 3.x and installing requirements.txt.

**When you build it**: Only when you need to add a new package to requirements.txt or upgrade the Python version.

**RPy2**: The "rpy2" package is actually built in the final Docker (umgear) image, due to some dependencies on R.

**The output**: This is currently built and pushed as adkinsrs/gear-python-base:YYYY-MM-DD

#### Dockerfile.r (The R Base)

This file is dedicated entirely to compiling R and running your Bioconductor scripts.

**When you build it**: Almost never. Only touch this if the team specifically requests a new version of Bioconductor or a brand-new R system library.

**The output**: This is currently built and pushed as adkinsrs/gear-r-base:YYYY-MM-DD

#### Dockerfile (The Final App)

This is your main daily-driver file. It starts with a clean Ubuntu image, uses COPY --from=... to pull in the pre-compiled folders from your registry, installs Apache, and copies over your Flask API and HTML/JS files.

**When you build it**: Every time you update the website, tweak the Apache configuration, or change a CGI script.  Anything gEAR-code related, basically.

**The output**: This builds in seconds and becomes your final production image.  This is pushed as adkinsrs/umgear:YYYY-MM-DD

## Starting the stack

If you used Method 1, in the `docker-compose.yml file` ensure the "build" step from the "web" service is commented out or deleted.

To start:
`docker compose up -d`
To stop:
`docker compose down -v`


Adding a service name (i.e. "web", "db") to the end of a command just performs this for that service.

### Get feature mapping files

Ask @adkinsrs for these (which will probably be in tar.gz format).  These hdf5 files should go in `<gear_root>/www/feature_mapping` so that orthology mapping will work.

## ProjectR

Currently the gear.ini.docker file is configured to send projectR jobs to the cloud, but to not use RabbitMQ.  This is to save me the hassle of managing an extra service.  What this means is that apache process will manage the jobs and job logging will write to /var/log/apache2/ssl_umgear_error.log in the container.

## "panel" service

If you do not plan on working on spatial panel stuff, feel free to comment out the "panel" service block from the docker-compose.yml file and skip this step.

You can pre-build the "panel_app" image using the Dockerfile from `<gear_root>/services/spatial/Dockerfile`

You can view logs with `docker compose logs panel`

## Getting datasets

Dataset files to use need to be housed on the host machine initially since they are not contained in the gEAR codebase and thus will not be available in the Dockerized version of gEAR out of the box

Generally, for development purposes, it is best to just have datasets for a couple of dataset collections, such as "Hearing (site default)" and "Ear (diverse variety)". The "Ear" collection does not exist on gEAR anymore but is a good representative collection of datasets.

The following dataset IDs are representative of both of the aforementioned dataset collections:

```text
7812a487-932b-32f7-2de7-33dd3155c849
c69485b2-6f8d-c60e-7337-e7ebad89b2c0
6fdd350c-4f82-07e2-3a39-408f105db16d
320ca057-0119-4f32-8397-7761ea084ed1
df726e89-b7ac-d798-83bf-2bd69d7f3b52
deb21a3b-677c-13e6-92cc-740fe8505e7c
bee735e5-d180-332c-7892-dd751dd76bb8
bad48d04-db27-26bc-2324-e88506f751fd
8779ce11-719d-58db-e626-c850e96a5379
09e5076e-754e-8738-30aa-5c7062ad9447
64485ca3-cf99-2993-99a3-54df3a09195c
cf8272cb-57fa-e841-0b50-9198e62fe2ff
2f4dc784-f581-6a43-0c51-0613b16c4930
1b12dde9-1762-7564-8fbd-1b07b750505f
a2dd9f06-5223-0779-8dfc-8dce7a3897e1
f7de7db2-b4cb-ebe3-7f1f-b278f46f1a7f
e34fa5c6-1083-cacb-eedf-23f59f2e005f
e78db16a-3927-348c-f5aa-4256330a7dff
f04627f7-6824-1e15-6a93-8de550c7b0a4
1531b46d-8435-f31a-3898-20abc2fc974a
cc90d264-ce9d-0b7b-1898-0e357c94b155
```

You can write some script to loop through these IDs and download them to your "datasets" directory. If your docker-compose.yml file mounting is configured correctly then you can do this outside of the container

`cd <gear_root>/www/datasets`
`wget https://umgear.org/datasets/<dataset_id>.h5ad .`

## MySQL

* See [MySQL Setup Notes](./docker_mysql.md) for more information.  Shaun has a dump file you can use (just ask @adkinsrs).  However do note that if there is a dataset or weighted gene list entry in the db and you do not have the physical file, then you may have errors.

## Running the docker container

The docker-compose.yml file is set up to mount the gEAR code as a volume, allowing you to make immediate edits to be reflected inside the container.  If making changes within `<gear_root>/www/api`, they will apply once you run `docker compose restart web`.

## Viewing logs

To view potential logs, run `docker compose logs` for all services or `docker compose logs <service>` for a single service.

If you want to view Apache logs, from server-side (Python) code, you can run `docker compose exec web tail -f /var/log/apache2/ssl_umgear_error.log` to view a running error log.  Sometimes it may be necessary to view "/var/log/apache2/error.log" instead.

## Issues and potential solutions

* I cannot log in
  * Clear your browser's cache. This can be quickly done with Ctrl-Shift-R (or Cmd-Shift-R on Mac)
  * Use the email address in the 'username' place.  I'm an idiot and always forget that.

* Running commands with executables give various errors related to library packages
  * If are you on a newer Mac OS system, and the default built image does not work, you may need to explicitly add `--platform=linux/amd64` to your build command options, to ensure the right libraries are being used.
  * If you change this, make sure the "web" service in the docker-compose.yml file has the `platform: linux/amd64` option added as well.
