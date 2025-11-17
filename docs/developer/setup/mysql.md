# MySQL configuration notes for gEAR instances

This document is intended to capture all the customizations to the mysql config files needed for a gEAR instance to operate properly.

## Installation

    sudo apt install mysql-server

### Setting a root password

This isn't entirely required, but I had to chase this issue down for some time.  On Ubuntu 18.04 and higher the root user is managed via the auth_socket plugin instead of mysql_native_password.  To fix this:

    mysql> ALTER USER 'root'@'localhost' IDENTIFIED WITH mysql_native_password BY 'whatever';
    Query OK, 0 rows affected (0.00 sec)

    mysql> flush privileges;
    Query OK, 0 rows affected (0.00 sec)

## Set up the initial database

For reference, the MySQL schema is housed [in this file](../../../create_schema.sql)

### New portals

```bash
     sudo mysql -u root -h localhost
```

You can change the name of the database as you like but need to make the corresponding change within the gear.ini file described below.

    mysql> create database gear_portal;
    mysql> use gear_portal;
    mysql> source /home/$USER/git/gEAR/create_schema.sql

### Existing portals setting up a duplicate instance

Assuming that you're running on GCP: use gcloud to grab the latest database dump

```bash
  cd /tmp
  gcloud compute scp gear-production:/home/jorvis/gear_portal.20191106.sql.gz .
  gunzip gear_portal.20191106.sql.gz
```

### The values here come from the gear.ini file in the repository root

  mysql> create database gear_portal;
  mysql> use gear_portal;
  mysql> source /tmp/gear_portal.20191106.sql

[ go run errands ... ]

  mysql> CREATE USER whoever@localhost IDENTIFIED BY 'whatever';
  mysql> GRANT select, insert, update, delete ON gear_portal.* TO whoever@localhost;
