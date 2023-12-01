# Set up the initial database

## Change the MySQL configuration to a wider character set by default

    $ vim /etc/mysql/mysql.conf.d/mysql.cnf

Add:

    [mysql]
    default-character-set=utf8mb4

Then:

    $ vim /etc/mysql/mysql.conf.d/mysqld.cnf

Add:

    [mysqld]
    character-set-server=utf8mb4
    collation-server=utf8mb4_unicode_ci

## New portals

    $ sudo mysql -u root -h localhost

You can change the name of the database as you like but need to make the corresponding change within the gear.ini file described below.

    mysql> create database gear_portal;
    mysql> use gear_portal;
    mysql> source /home/$USER/git/gEAR/create_schema.sql

## Existing portals setting up a duplicate instance:

Assuming that you're running on GCP: use gcloud to grab the latest database dump

$ cd /tmp
$ gcloud compute scp gear-production:/home/jorvis/gear_portal.20191106.sql.gz .
$ gunzip gear_portal.20191106.sql.gz

# the values here come from the gear.ini file in the repository root
mysql> create database gear_portal;
mysql> use gear_portal;
mysql> source /tmp/gear_portal.20191106.sql

[ go run errands ... ]

mysql> CREATE USER whoever@localhost IDENTIFIED BY 'whatever';
mysql> GRANT select, insert, update, delete ON gear_portal.* TO whoever@localhost;
