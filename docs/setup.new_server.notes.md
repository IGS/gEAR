# Info for setting up or upgrading a new server

## These operations were last performed on Ubuntu 22.04 LTS

Instances of a gEAR Portal are most often run within a cloud instance, where you can choose your own operating system and resources.  On Google Cloud for a starter instance I chose an e2-standard-2 (2 vCPUs and 48GB RAM) with 300GB of solid state disk space.  You'll definitely want to increase the CPU as you gain more simultaneous users and RAM depending on your dataset sizes.  Once you create and start the instance:

```bash
    sudo apt update
    sudo apt upgrade
```

Reboot if there are kernel updates (or just to be safe if you don't know.)

```bash
    cd && mkdir git
    sudo apt install git
    cd git
    git clone https://github.com/IGS/gEAR.git
```

### MYSQL

    sudo apt install mysql-server

Follow instructions in our setup.mysql.md document

### R

Not necessary if you want projectR to run on a Google Cloud Run service (configurable in gear.ini)

Please consult `setup.r_rpy2.md` for packages to install in order to install R and requisite R packages

### RabbitMQ

Not necessary if you want projectR to run in the Apache environment or do not want to setup the RabbitMQ messaging service (configurable in gear.ini)

Follow instructions in setup.rabbit_mq.md document

### Python

Follow instructions in setup.python.md document

### APACHE

    sudo apt install apache2 apache2-dev

Follow instructions in setup.apache.md document

### gEAR portal

```bash
cd /var
sudo rm -rf www && sudo ln -s ~/git/gEAR/www
```

### Data transfer

If moving to a new instance and you need to transfer data you need
to shift the contents of the following directories:

$HOME/git/gEAR/www/datasets
$HOME/git/gEAR/www/analyses

### Permissions

Permissions need to be writeable for your apache user.

```bash
cd /var/www
chmod 777 datasets analyses/* uploads/files/
```

### (For devel servers) Github Actions self-hosted runner

The idea here is that when commiting code, we can use Github Actions to pull the code down on to the server to do things like automated testing.

(to be continued in https://github.com/IGS/gEAR/issues/632)

### Undocumented currently

- Setting up mail server so forgot password works
- What to do for helpdesk config?
