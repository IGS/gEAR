# MySQL Docker setup notes

To be performed after performing "docker-compose up -d"

## Get root password

Do the following to get the root password:
`docker-compose logs db | grep "GENERATED"`

This will get you the randomly generated root password.  Note that the password will be different every time the docker-compose stack is brought up.

## Pull some backup data into the mysql container (first-time only)

This step only needs to be performed the first time you create the MySQL container on your computer with "docker-compose".  If you have a "mysql" directory already present within this "docker" directory, then this step is most likely not necessary.  This is because in the docker-compose.yml file, the "mysql" service mounts a volume from "/var/lib/mysql" in the container to "./mysql" (assuming the `docker-compsoe up -d` command was run in the same directory as this file).

To pull backup file to your host machine, do:
`gcloud compute scp gear-devel:/home/jorvis/gear_portal.20191106.sql.gz .`

This file, when gunzipped will be about 1.6 Gb.  After this, we copy the file into the Docker "mysql" container.

1. Uncompress the file with `gunzip gear_portal.20191106.sql.gz`
2. Figure out the container ID with `docker ps`.  It should be the first field in the result row where the "IMAGE" is "mysql:5.7"
3. Run `docker cp ./gear_portal.20191106.sql <container_id>:/tmp/gear_portal.20191106.sql

## Set up the mysql database (first-time only)

1. Do `docker-compose exec mysql /bin/bash` to get into the docker instance.  Next do `mysql -uroot -p<GENERATED_ROOT_PASSWORD> gear\_portal` to log into mysql as root.  Note that the "GENERATED_ROOT_PASSWORD" was obtained from the "Get root password" section, and that there is no space between the "-p" and the password.
2. In the mysql client, run `source /tmp/gear\_portal.20191106.sql` to load the SQL dump backup file.
3. After that finishes run `grant select, insert, update, delete on gear\_portal.\* to gear@localhost identified by 'gearadmin';` to set up privileges for the "gear" user.

## Issues

### Cannot log in

Make sure the gear.ini file in the gEAR root directory has the host entry as "db" instead of "localhost".