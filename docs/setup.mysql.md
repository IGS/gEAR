# use gcloud to grab the latest database dump
$ cd /tmp
$ gcloud compute scp gear-production:/home/jorvis/gear_portal.20191106.sql.gz .
$ gunzip gear_portal.20191106.sql.gz

# the values here come from the gear.ini file in the repository root
mysql> create database gear_portal;
mysql> use gear_portal;
mysql> source /tmp/gear_portal.20191106.sql

[ go run errands ... ]

mysql> CREATE USER gear@localhost IDENTIFIED BY 'gearadmin';
mysql> GRANT select, insert, update, delete ON gear_portal.* TO gear@localhost;
