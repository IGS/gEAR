# Notes about the Docker setup

* To build run `docker build -t umgear:devel --build-arg git_user=<user> --build-arg git_pass=<pass> .`
  * If you tag it under a new image, ensure it is reflected in the docker-compose.yml file
* Changed some hardcoded paths that pointed to "/home/jorvis/git/gEAR" to reflect the relative path of the executed file looking for that path.  This should make it compatible with both Docker with devel/production versions of the codebase.
* Later on I would love to do a multi-stage build of the Dockerfile so that the "-dev" ubuntu packages are not in the final container... only compiled Ubuntu and Python packages we need.

## Annotations and Datasets

* Annotation and dataset files to use need to be housed on the host machine initially since they are not contained in the gEAR codebase and thus will not be available in the Dockerized version of gEAR out of the box
* Right now the only annotation I am including is Mouse - release 94.  This is because each annotation directory is about 1Gb each and would overrun the storage on my computer and blow up the size of the 'dataset' Docker image if all were added.  We only need one annotation for sandboxing and development though.
* I am also using the Hertzano/Ament P2 mouse cochlea dataset for testing (as recommended by Brian Herb)

### Loading annotation data

* After launching the Docker containers via docker-compose, the annotation data may not be loaded into the MySQL database.  To do this (substitute species and release number as needed):

1. Run `docker-compose exec web /bin/bash`
2. `cd /opt/gEAR-dockerized/annotations/mouse`
3. Run `../../bin/load_genbank_annotations.py -i ./release-94 -id 1 -r 94` which will load the annotation file

## Issues and potential solutions

* I cannot log in
  * Clear your browser's cache
  * Use the email address in the 'username' place.  I'm an idiot and always forget that.