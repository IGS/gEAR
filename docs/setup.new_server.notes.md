# Info for setting up or upgrading a new server

## These operations were last performed on Ubuntu 20.04 LTS

Instances of a gEAR Portal are most often run within a cloud instance, where you can choose your own operating system and resources.  On Google Cloud for a starter instance I chose an e2-standard-2 (2 vCPUs and 48GB RAM) with 300GB of solid state disk space.  You'll definitely want to increase the CPU as you gain more simultaneous users and RAM depending on your dataset sizes.  Once you create and start the instance:

    $ sudo apt update
    $ sudo apt upgrade

Reboot if there are kernel updates (or just to be safe if you don't know.)

    $ cd && mkdir git
    $ sudo apt install git
    $ git clone https://github.com/jorvis/gEAR.git

### Updating Ubuntu version to 22.04 LTS

To check current version, perform `lsb_release -a`

`sudo apt install update-manager-core` (which should already be installed)
`sudo apt update && sudo apt dist-upgrade`
`sudo do-release-upgrade`

You may be prompted to perform a reboot at this point in order to do the upgrade, which is done with `sudo reboot`.  This will kick you out of the VM.  Just ssh back in and do `sudo do-release-upgrade`.

Follow the prompts and let the upgrade do its thing. There is a note that it can take several hours, so keep that in mind.  The upgrade will also prompt for another restart of the server.

Do another `lsb_release -a` to confirm the version upgrade (should be 22.04 Jammy)

`sudo apt update`

At this point, we could remove the Python2 packages using `sudo apt autoremove`. From my experience Python2-related packages were the only things to go.

### MYSQL

    $ sudo apt install mysql-server

Follow instructions in our setup.mysql.md document

### R

Not necessary if you want projectR to run on a Google Cloud Run service (configurable in gear.ini)

`sudo apt install r-base`

Please consult `setup_notes_r_rpy2.md` for packages to install in order to install requisite R packages

### RabbitMQ

Not necessary if you want projectR to run in the Apache environment or do not want to setup the RabbitMQ messaging service (configurable in gear.ini)

Follow instructions in setup_rabbit_mq.md document

### Python

Follow instructions in setup.python.md document

### APACHE

    $ sudo apt install apache2 apache2-dev

Follow instructions in setup.apache.md document

### gEAR portal

```bash
$ cd ~/git
$ git clone https://github.com/jorvis/gEAR.git
$ cd /var
$ sudo rm -rf www && sudo ln -s ~/git/gEAR/www
```

### Data transfer

If moving to a new instance and you need to transfer data you need
to shift the contents of the following directories:

$HOME/git/gEAR/www/datasets
$HOME/git/gEAR/www/analyses

### Permissions

Permissions need to be writeable for your apache user.

```bash
$ cd /var/www
$ chmod 777 datasets analyses/* uploads/files/
```
