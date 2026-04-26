# Info for setting up or upgrading a new server

## These operations were last performed on Ubuntu 22.04 LTS

Instances of a gEAR Portal are most often run within a cloud instance, where you can choose your own operating system and resources.  On Google Cloud for a starter instance I chose an e2-standard-2 (2 vCPUs and 48GB RAM) with 300GB of solid state disk space.  You'll definitely want to increase the CPU as you gain more simultaneous users and RAM depending on your dataset sizes.  Once you create and start the instance:

NOTE: Extra apt install commands are present in the R setup and in the python setup.

```bash
    sudo apt update
    sudo apt upgrade
    sudo apt install build-essential \  # installs gcc, g++, make
        curl \
        rsync \
        vim \   # Shaun Adkins was here
        unzip \
        wget \
        ca-certificates \
        fonts-roboto \    # Matplotlib font
        ttf-mscorefonts-installer \    # Matplotlib font alternative
        fontconfig \
    sudo fc-cache -f -v    # cache the fonts
```

Reboot if there are kernel updates (or just to be safe if you don't know.)

```bash
    cd && mkdir git
    sudo apt install git
    cd git
    git clone https://github.com/IGS/gEAR.git
```

### MYSQL

`sudo apt install mysql-server`

Follow instructions in [the MySQL setup](./mysql.md) document

### R

Not necessary if you want projectR to run on a Google Cloud Run service (configurable in gear.ini)

Please consult [the R setup](./r_rpy2.md) for packages to install in order to install R and requisite R packages

### RabbitMQ

Not necessary if you want projectR to run in the Apache environment or do not want to setup the RabbitMQ messaging service (configurable in gear.ini)

Follow instructions in [the RabbitMQ setup](./rabbitmq.md) document

### Python

Follow instructions in [the Python setup](./python.md) document

### APACHE

`sudo apt install apache2 apache2-dev`

Follow instructions in [the Apache setup](./apache.md) document

### Sass

Used for changing theme colors between portal flavors.
NOTE: Ruby Sass is end-of-life and they recommend switching to Dart Sass

`sudo apt install ruby-sass`

### gEAR portal

```bash
cd ~jorvis/git
git clone https://github.com/jorvis/gEAR.git
cd /var
sudo rm -rf www && sudo ln -s ~jorvis/git/gEAR/www
```

### Systemd Services

More information about these services can be found at `gEAR/systemd/README.md`

```bash
cd ~jorvis/git/gEAR/systemd
sudo cp *target /etc/systemd/system/
sudo cp *service /etc/systemd/system/

echo "**IMPORTANT**: For the *service files, correct the <gEAR_path> to point to the gEAR root on this server"

# Start the services
cd /etc/systemd/system

sudo systemctl enable projectr-consumer.target gosling-upload-consumer.target
sudo systemctl start projectr-consumer.target gosling-upload-consumer.target

sudo systemctl enable spatial-panel.service
sudo systemctl start spatial-panel.service

# <Add more services here>
```

### Data transfer

If moving to a new instance and you need to transfer data you need
to shift the contents of the following directories:

$HOME/git/gEAR/www/datasets
$HOME/git/gEAR/www/analyses

### Permissions

Permissions need to be writeable for the apache user (www-data).

```bash
cd /var/www
chmod 777 datasets datasets/spatial analyses/* carts/ projections/ uploads/files/ img/dataset_previews/ cache/spatial_panel/
```

### Undocumented currently

- Setting up mail server so forgot password works
- What to do for helpdesk config?
