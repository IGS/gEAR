## These operations were last performed on Ubuntu 20.04 LTS

Instances of a gEAR Portal are most often run within a cloud instance, where you can choose your own operating system and resources.  On Google Cloud for a starter instance I chose an e2-standard-2 (2 vCPUs and 8GB RAM) with 300GB of solid state disk space.  You'll definitely want to increase this as you gain more simultaneous users.  Once you create and start the instance:

    $ sudo apt update
    $ sudo apt upgrade

Reboot if there are kernel updates (or just to be safe if you don't know.)

    $ cd && mkdir git
    $ sudo apt install git
    $ git clone https://github.com/jorvis/gEAR.git

### MYSQL

    $ sudo apt install mysql-server

Follow instructions in our setup.mysql.md document

### Python

Follow instructions in setup.python.md document

### APACHE

    $ sudo apt install apache2 apache2-dev

Follow instructions in setup.apache.md document

### gEAR portal

    $ cd ~/git
    $ git clone https://github.com/jorvis/gEAR.git
    $ cd /var
    $ sudo rm -rf www && sudo ln -s ~/git/gEAR/www

### Data transfer
#
#  If moving to a new instance and you need to transfer data you need
#  to shift the contents of the following directories:

$HOME/git/gEAR/www/datasets
$HOME/git/gEAR/www/analyses

### Permissions

Permissions need to be writeable for your apache user.

    $ cd /var/www
    $ chmod 777 datasets analyses/* uploads/files/

