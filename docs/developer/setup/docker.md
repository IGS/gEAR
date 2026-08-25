# Notes about the Docker setup

## Before building

* [Orbstack](https://orbstack.dev/) is recommended over Docker Desktop, but certainly optional. It is faster, lighter, and still uses the same Docker command-line tools.  Docker Desktop gets the job done just fine as well.
  * Unfortunately Orbstack does not have Windows support.
* Perform a `git clone` on the gEAR repository if you have not.  This is so you can hot-load code changes into your container and still be able to commit it.
* From the gEAR root, `cd docker`. All commands assume you are in this directory
* `cp docker-compose.yml.template docker-compose.yml`
  * Alternatively ask @adkinsrs for a docker-compose.yml file as it will be filled in. Otherwise fill in any values wrapped in brackets
  * @adkinsrs's file is hard-coded to his paths so be sure to change those.
* `mkdir mysql` - Necessary to mount the MySQL data directory for preserving

## Files from @adkinsrs

There are various .template files hanging around with some missing paths.  @adkinsrs has live versions of the files that he can give you.

* gear.ini - Goes in <gear_root> directory (parent directory of "docker")
* docker-compose.yml - Goes in "docker" directory. Ensure the volume mount paths reflect your own gEAR repo location.
* feature_mapping.tar.gz - Extract this directory into `<gear_root>/www/`. Contains orthology mapping files in hdf5 format.
* MySQL dump file.  Please consult [MySQL Setup Notes](./docker_mysql.md) to load this file after you bring the docker compose stack up
  * Do note that if there may be weighted gene list entries in the db and you do not have the physical file, then you may have errors.

### Getting datasets

Dataset files to use need to be housed on the host machine initially since they are not contained in the gEAR codebase and thus will not be available in the Dockerized version of gEAR out of the box

Generally, for development purposes, it is best to just have datasets for a couple of dataset collections, such as "Hearing (site default)".

The following dataset IDs are representative of both of the aforementioned dataset collections:

```text
7812a487-932b-32f7-2de7-33dd3155c849
c69485b2-6f8d-c60e-7337-e7ebad89b2c0
6fdd350c-4f82-07e2-3a39-408f105db16d
8779ce11-719d-58db-e626-c850e96a5379
cf8272cb-57fa-e841-0b50-9198e62fe2ff
deb21a3b-677c-13e6-92cc-740fe8505e7c
64485ca3-cf99-2993-99a3-54df3a09195c
320ca057-0119-4f32-8397-7761ea084ed1
df726e89-b7ac-d798-83bf-2bd69d7f3b52
09e5076e-754e-8738-30aa-5c7062ad9447
2f4dc784-f581-6a43-0c51-0613b16c4930
bee735e5-d180-332c-7892-dd751dd76bb8
bad48d04-db27-26bc-2324-e88506f751fd
```

You can write some script to loop through these IDs and download them to your "datasets" directory. If your docker-compose.yml file mounting is configured correctly then you can do this outside of the container

`cd <gear_root>/www/datasets`
`wget https://umgear.org/datasets/<dataset_id>.h5ad .`
`cd <gear_root>/www/datasets_uploaded`
`wget https://umgear.org/datasets_uploaded/<dataset_id>.svg .`

## Acquiring the gEAR images

The easiest way to get the gEAR stack is to `cd docker` (if you have not done so) and run `docker compose up -d`, which will pull the images off of Docker Hub.

Running `docker compose up -d` will also perform a check to see if the image (typically the "latest" tag) has been updated in the Docker Hub registry and will pull that down if it needs to.

After all containers are up, if you have not done so set up MySQL [MySQL Setup Notes](./docker_mysql.md)

## Starting the stack

To start:
`docker compose up -d`

Perform `docker compose ps` to see the status of running containers.

To stop:
`docker compose down -v`

You can also pass the name of a specific service, like "web" to any docker compose command to only do this just for that service.  An example would be `docker compose up -d web`

The docker-compose.yml file is set up to mount the gEAR code as a volume, allowing you to make immediate edits to be reflected inside the container.  If making changes within `<gear_root>/www/api`, they will apply once you run `docker compose restart web`.

Adding a service name (i.e. "web", "db") to the end of a command just performs this for that service.

### Viewing logs

To view potential logs, run `docker compose logs` for all services or `docker compose logs <service>` for a single service.

If you want to view Apache logs, from server-side (Python) code, you can run `docker compose exec web tail -f /var/log/apache2/ssl_umgear_error.log` to view a running error log.  Sometimes it may be necessary to view "/var/log/apache2/error.log" instead.

## Other services

### ProjectR

Currently the gear.ini.docker file is configured to send projectR jobs to the cloud, but to not use RabbitMQ.  This is to save me the hassle of managing an extra service.  What this means is that apache process will manage the jobs and job logging will write to /var/log/apache2/ssl_umgear_error.log in the container.

### "panel" service

If you do not plan on working on spatial panel stuff, feel free to comment out the "panel" service block from the docker-compose.yml file and skip this step.

You can pre-build the "panel_app" image using the Dockerfile from `<gear_root>/services/spatial/Dockerfile`

You can view logs with `docker compose logs panel`

## Building the images manually (Advanced)

These steps will both build all images and push to the adkinsrs Docker Hub repository.  You will not have access to this repo (unless you are Shaun), so you just do the `docker compose up -d` method, or rewrite `Dockerfile`, `docker-bake.hci`, and `docker-compose.yml` to point to your own space.

1. `cd docker`
2. Build intermediate images with `docker buildx bake --allow=fs.read=.. intermediate`
3. Build the docker compose stack images with `docker buildx bake --allow=fs.read=.. default`

### The three "umgear" Dockerfiles

Information about the three Dockerfiles found in `<gear_root>/docker`

#### Dockerfile.python (The Python Base)

This file is dedicated entirely to compiling Python 3.x and installing requirements.txt.

**When you build it**: Only when you need to add a new package to requirements.txt or upgrade the Python version.

**RPy2**: The "rpy2" package is actually built in the final Docker (umgear) image, due to some dependencies on R.

**The output**: This is currently built and pushed as adkinsrs/gear-python-base:YYYY-MM-DD and also tagged with the "latest" tag.

#### Dockerfile.r (The R Base)

This file is dedicated entirely to compiling R and running your Bioconductor scripts.

**When you build it**: Almost never. Only touch this if the team specifically requests a new version of Bioconductor or a brand-new R system library.

**The output**: This is currently built and pushed as adkinsrs/gear-r-base:YYYY-MM-DD and also tagged with the "latest" tag.

#### Dockerfile (The Final App)

This is your main daily-driver file. It starts with a clean Ubuntu image, uses COPY --from=... to pull in the pre-compiled folders from your registry, installs Apache, and copies over your Flask API and HTML/JS files.

Currently the inherited R and Python images are set to use the "latest" tag of a locally built image, as most of the time we want the most up-to-date version. If for some reason you need an earlier version, edit the Dockerfile to use one of the existing `adkinsrs/<image>:YYYY-MM-DD` tags stored in Docker Hub.

**When you build it**: Every time you update the website, tweak the Apache configuration, or change a CGI script.  Anything gEAR-code related, basically.

**The output**: This builds in seconds and becomes your final production image.  This is pushed as adkinsrs/umgear:YYYY-MM-DD and also tagged with the "latest" tag.

## Issues and potential solutions

* I cannot log in
  * Clear your browser's cache. This can be quickly done with Ctrl-Shift-R (or Cmd-Shift-R on Mac)
  * Use the email address in the 'username' place.  I'm an idiot and always forget that.

* Running commands with executables give various errors related to library packages
  * If are you on a newer Mac OS system, and the default built image does not work, you may need to explicitly add `--platform=linux/amd64` to your build command options, to ensure the right libraries are being used.
  * If you change this, make sure the "web" service in the docker-compose.yml file has the `platform: linux/amd64` option added as well.
