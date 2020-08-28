# Intro

This document explains the process to setup Epiviz with gEAR/NEMO


# Epiviz API

the API is setup in `www/epiviz-api` folder


## Installing dependencies

follow the instructions in the README to install the dependencies

```
# create a virtual env
virtualenv env --python=/opt/bin/python3

# install dependencies
pip install -r requirements.txt
```

## Supervisor setup

Epiviz uses Sanic, an asynchronous web server not supported by WSGI (needs ASGI, but apache doesn't support it directly yet). So we use gunicorn + supervisor to set things up. 

- Install these packages

    `sudo apt install gunicorn supervisor`

- Add the gunicorn config file so that supervisor manages this process (create this file `/etc/supervisor/conf.d/epiviz.conf`)

    ```
        [program:gunicorn]
        directory=/var/www/epiviz-api
        environment=PYTHONPATH=/var/www/epiviz-api/bin/python
        command=/var/www/epiviz-api/env/bin/gunicorn epiviz:app --log-level debug --bind 0.0.0.0:8000 --worker-class sanic.worker.GunicornWorker
        autostart=true
        autorestart=true
        stderr_logfile=/var/log/gunicorn/gunicorn.err.log
        stdout_logfile=/var/log/gunicorn/gunicorn.out.log
    ```


- Let supervisor run and update this process

    ```
        supervisorctl reread
        supervisorctl update

        service supervisor restart
    ```


## Add Proxy Endpoints to Apache Configuration

Add this to the end of `/etc/apache2/sites-available/umgear-ssl.conf` before the `</VirtualHost>` closing tag

    ProxyPreserveHost On
    <Location "/epiviz">
    ProxyPass "http://127.0.0.1:8000/"
    ProxyPassReverse "http://127.0.0.1:8000/"
    </Location>

    <Location "/epivizfile">
    ProxyPass "http://127.0.0.1:8000/file/"
    ProxyPassReverse "http://127.0.0.1:8000/file/"
    </Location>

    <Location "/trackhub">
    ProxyPass "http://127.0.0.1:8000/trackhub/"
    ProxyPassReverse "http://127.0.0.1:8000/trackhub/"
    </Location>

    <Directory "/datasets_epigenetic">
    Options Indexes FollowSymLinks
    </Directory>

Restart Apache

    `service apache2 restart`


# File upload instructions

All epigenetic data uploaded to the gEAR interface are stored in `datasets_epigenetic` folder

The current setup is to use a separate drive to store these files. IF thats the case, 

- Follow instructions on google cloud to attach a new disk to VM available @ https://cloud.google.com/compute/docs/disks/add-persistent-disk
- move the datasets_epigenetic directory to the new disk & add a simlink to the drive
    
    if the mount point is `/mnt/disks/epiviz-drive`

    ```
        mv datasets_epigenetic /mnt/disks/epiviz-drive
        ln -s /mnt/disks/epiviz-drive datasets_epigenetic
    ```

- Apache requires this folder to have full permissions

    `chmod -R 777 datasets_epigenetic`

- Mount the directory on instance reboots

    `sudo mount -o discard,defaults /dev/sdb /mnt/disks/epiviz-drive`


# Add epiviz dataset table to the gEAR/NEMO DB

(also added to `create_scheme.sql`)

    ```
        CREATE TABLE `dataset_epiviz` (
        `id` varchar(50) NOT NULL,
        `owner_id` int(11) NOT NULL,
        `annotation` text,
        `type` varchar(10) NOT NULL,
        `url` text NOT NULL,
        `title` text,
        `is_public` tinyint(4) DEFAULT '0',
        `description` text,
        `share_id` varchar(50) NOT NULL,
        `organism` varchar(100) DEFAULT NULL
        ) ENGINE=InnoDB;
    ```