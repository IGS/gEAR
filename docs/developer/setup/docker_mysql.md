# MySQL Docker setup notes

To be performed after performing "docker-compose up -d"

## Get root password

Root password is set in the docker-compose.yml file.

## Pull some backup data into the mysql container (first-time only)

This step only needs to be performed the first time you create the MySQL container on your computer with "docker-compose".  If you have a "mysql" directory already present within this "docker" directory, then this step is most likely not necessary.  This is because in the docker-compose.yml file, the "mysql" service mounts a volume from "/var/lib/mysql" in the container to "./mysql" (assuming the `docker-compose up -d` command was run in the same directory as this file).

To pull backup file to your host machine, do:
`gcloud compute scp <server>:<db_dump.sql> .`

Alternative ask @adkinsrs for a SQL dump file.

This file, when gunzipped will be about 1.6 Gb.  After this, we copy the file into the Docker "mysql" container.

1. Uncompress the file with `gunzip` if you need to.
2. Figure out the container ID with `docker ps`.  It should be the first field in the result row where the "IMAGE" is "mysql:8.0"
3. Run `docker cp ./<db_dump.sql> <container_id>:/tmp/<db_dump.sql>`

NOTE: Change the SQL filename to whatever database dump you are using.

## Set up the mysql database

### Dump file

1. Do `docker-compose exec db /bin/bash` to get into the docker instance.  Next do `mysql -uroot -p<ROOT_PASSWORD>` to log into mysql as root.  Note that the "GENERATED_ROOT_PASSWORD" was obtained from the "Get root password" section, and that there is no space between the "-p" and the password.
1. In the mysql client, run `source <db_dump.sql>` to load the SQL dump backup file.

### No dump file (fresh container only)

1. Do `docker-compose exec db /bin/bash` to get into the docker instance.  Next do `mysql -uroot -p<ROOT_PASSWORD>` to log into mysql as root.  Note that the "GENERATED_ROOT_PASSWORD" was obtained from the "Get root password" section, and that there is no space between the "-p" and the password.
2. Run the following (in the mysql client) to setup the initial MySQL tables.
    1. `create database gear_portal;`
    2. `use gear_portal;`
    3. `source <gear_root> create_schema.sql`
4. After that finishes run the following to ensure the gEAR user can do database operations in gEAR
    1. `GRANT USAGE ON *.* TO 'gear'@'%';`
    2. `GRANT SELECT, INSERT, UPDATE, DELETE ON gear_portal.* TO 'gear'@'%';`

## Issues

### Cannot log in

Make sure the gear.ini file in the gEAR root directory has the host entry as "db" instead of "localhost".

### MySQL container will not start

If you check `docker compose logs db` and it says something about `chown: cannot dereference '/var/lib/mysql/mysql.sock': No such file or directory`, just delete the `./mysql/mysql.sock` file in this directory, then do `docker compose down -v; docker compose up -d`. It should start up properly with a new socket file