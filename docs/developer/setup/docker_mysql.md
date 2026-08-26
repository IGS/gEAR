# MySQL Docker setup notes

To be performed after performing `docker-compose up -d`. The service name is called "db".

For reference, the MySQL schema is housed [in this file](../../../create_schema.sql)

## Get root password

Root password is set in the docker-compose.yml file under `MYSQL_ROOT_PASSWORD`

## Create a MySQL dump file

This ideally will only be once, or none if you are given a dump file by @adkinsrs

Because our MySQL servers store a lot of dataset and gene metadata, the dump file can be large.  Fortunately I made a script to reduce the size as much as possible.

1. `docker compose exec web /bin/bash`
2. `cd /opt/gEAR/bin`
3. `/opt/bin/python3 create_test_mysql_dump.py --layout-ids 0 --dump-file /tmp/mini-gear.sql --user root --pass <ROOT_PASSWORD>`
    1. These layout_ids corresponsd to the "Hearing (site default)"
    2. This is also set to only pull the genes from the "gene" table that are the most recent ensembl release version for each organism ID
    3. The user and pass in the command line are specific to the docker mysql service. Using the default user and pass will not work as the default user does not have PROCESS privileges to dump.
4. `docker compose cp web:/tmp/mini-gear.sql ./gear-mini.sql` to copy off the docker container

## Pull some backup data into the mysql container (first-time only)

This step only needs to be performed the first time you create the MySQL container on your computer with "docker-compose".  If you have a "mysql" directory already present within this "docker" directory, then this step is most likely not necessary.  This is because in the docker-compose.yml file, the "mysql" service mounts a volume from "/var/lib/mysql" in the container to "./mysql" (assuming the `docker-compose up -d` command was run in the same directory as this file).

Run `docker compose cp ./<db_dump.sql> db:/tmp/<db_dump.sql>`. This sets you up for Method 1

If you do not have a dump file, just use the `<gear_root>/create_schema.sql` file and Method 2 below.

NOTE: Change the SQL filename to whatever database dump you are using.

## Set up the mysql database

### Method 1: From a dump file

1. Do `docker-compose exec db /bin/bash` to get into the docker instance.  Next do `mysql -uroot -p<ROOT_PASSWORD>` to log into mysql as root.  Note that the "GENERATED_ROOT_PASSWORD" was obtained from the "Get root password" section, and that there is no space between the "-p" and the password.
2. Run the following (in the mysql client) to setup the initial MySQL tables.
    1. `create database gear_portal;`
    2. `use gear_portal;`
    3. `source /tmp/<db_dump.sql>` to load the SQL dump backup file.  This will overwrite the database currently being used.
3. After that finishes run the following to ensure the gEAR user can do database operations in gEAR
    1. `GRANT USAGE ON *.* TO 'gear'@'%';`
    2. `GRANT SELECT, INSERT, UPDATE, DELETE ON gear_portal.* TO 'gear'@'%';`

### Method 2; No dump file (fresh container only)

1. Do `docker-compose exec db /bin/bash` to get into the docker instance.  Next do `mysql -uroot -p<ROOT_PASSWORD>` to log into mysql as root.  Note the root password from the "Get root password" section, and that there is no space between the "-p" and the password.
2. Run the following (in the mysql client) to setup the initial MySQL tables.
    1. `create database gear_portal;`
    2. `use gear_portal;`
    3. `source /tmp/create_schema.sql`
3. After that finishes run the following to ensure the gEAR user can do database operations in gEAR
    1. `GRANT USAGE ON *.* TO 'gear'@'%';`
    2. `GRANT SELECT, INSERT, UPDATE, DELETE ON gear_portal.* TO 'gear'@'%';`

## Issues

### Cannot log in

Make sure the gear.ini file in the gEAR root directory has the host entry as "db" instead of "localhost".

### MySQL container will not start

If you check `docker compose logs db` and it says something about `chown: cannot dereference '/var/lib/mysql/mysql.sock': No such file or directory`, just delete the `./mysql/mysql.sock` file in this directory, then do `docker compose down -v; docker compose up -d`. It should start up properly with a new socket file
