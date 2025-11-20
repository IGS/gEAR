# Notes about the Docker setup

## Before building

* From the gEAR root, `cd docker`
* `cp gear.ini.docker.template gear.ini.docker`
  * Alternatively ask @adkinsrs for a gear.ini.docker file as it will be filled in. Otherwise fill in any values wrapped in brackets
* `cp docker-compose.yml.template docker-compose.yml`
  * Alternatively ask @adkinsrs for a docker-compose.yml file as it will be filled in. Otherwise fill in any values wrapped in brackets
  * @adkinsrs's file is hard-coded to his paths so be sure to change those.

## Acquiring the gEAR image

There are two options here.  The first method is significantly quicker, but there is a chance it may not be updated based on the latest Dockerfile instructions (if @adkinsrs forgets to push the latest image up).  The first method may also not be possible if you are restricted by Google Cloud authentication for the gEAR project area.

### Method 1: Pull image

* Authenticate into Google Cloud
  * `gcloud auth configure-docker us-east1-docker.pkg.dev`
* Pull the image
  * `docker pull us-east1-docker.pkg.dev/gear-154704/gear-image/umgear`
* Tag the image
  * `docker tag us-east1-docker.pkg.dev/gear-154704/gear-image/umgear umgear:main`

### Method 2: Build image

* Ensure you are in the "devel" branch of gEAR before building (`git checkout devel`)
* To build run `docker build -t umgear:main .`
  * If you tag it under a new image, ensure it is reflected in the docker-compose.yml file
* The build can take a while, particularly in the Bioconductor installation steps. Fortunately completed steps are cachable.

In the build, the "gear.ini.docker" file will end up copied to "gear.ini" in the "/opt/gEAR" directory for the docker instance. However, if are using docker-compose and the gEAR directory is mounted into the "web" service, this can be overriden to a gear.ini from outside.  If you do not have a "gear.ini" file (only gear.ini.template), then ask @adkinsrs for one.

## Starting the stack

To start:
`docker compose up -d`
To stop:
`docker compose down -v`

Adding a service name (i.e. "web", "db") to the end of a command just performs this for that service.

## ProjectR

Currently the gear.ini.docker file is configured to send projectR jobs to the cloud, but to not use RabbitMQ.  This is to save me the hassle of managing an extra service.  What this means is that apache process will manage the jobs and job logging will write to /var/log/apache2/ssl_umgear_error.log in the container.

## "panel" service

You can pre-build the "panel_app" image using the Dockerfile from `<gear_root>/services/spatial/Dockerfile`

You can view logs with `docker compose logs panel`

## Annotations and Datasets

* Annotation and dataset files to use need to be housed on the host machine initially since they are not contained in the gEAR codebase and thus will not be available in the Dockerized version of gEAR out of the box
* Right now the only annotation I am including is Mouse - release 94.  This is because each annotation directory is about 1Gb each and would overrun the storage on my computer and blow up the size of the 'dataset' Docker image if all were added.  We only need one annotation for sandboxing and development though.
* I am also using the Hertzano/Ament P2 mouse cochlea dataset for testing (as recommended by Brian Herb)

### Loading annotation data

* After launching the Docker containers via docker-compose, the annotation data may not be loaded into the MySQL database.  To do this (substitute species and release number as needed):

1. Run `docker compose exec web /bin/bash`
2. `cd /opt/gEAR/annotations/mouse`
3. Run `../../bin/load_genbank_annotations.py -i ./release-94 -id 1 -r 94` which will load the annotation file

## Getting datasets

Generally, for development purposes, it is best to just have datasets for a couple of dataset collections, such as "Hearing (site default)" and "Ear (diverse variety)"

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
```

You can write some script to loop through these IDs and download them to your "datasets" directory. If your docker-compose.yml file mounting is configured correctly then you can do this outside of the container

`cd <gear_root>/www/datasets`
`wget https://umgear.org/datasets/<dataset_id>.h5ad .`

## MySQL

* See [MySQL Setup Notes](./docker_mysql.md) for more information.  Shaun has a dump file you can use (just ask @adkinsrs).  However do note that if there is a dataset or weighted gene list entry in the db and you do not have the physical file, then you may have errors.

## Running the docker container

The docker-compose.yml file is set up to mount the gEAR code as a volume, allowing you to make immediate edits to be reflected inside the container.  If making changes within `<gear_root>/www/api`, they will apply once you run `docker compose restart web`.

At any time you can run `docker compose exec web tail -f /var/log/apache2/ssl_umgear_error.log` to view a running error log.

## Spatial panel

```
cd <gear_root>/spatial
docker build -t panel_app .
```

This will enable the spatial panel dashboard to be used in the docker-compose.yml stack

## Issues and potential solutions

* I cannot log in
  * Clear your browser's cache
  * Use the email address in the 'username' place.  I'm an idiot and always forget that.
