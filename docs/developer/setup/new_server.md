# Info for setting up or upgrading a new server

## These operations were last performed on Ubuntu 22.04 LTS

Instances of a gEAR Portal are most often run within a cloud instance, where you can choose your own operating system and resources.  On Google Cloud for a starter instance I chose an e2-standard-2 (2 vCPUs and 48GB RAM) with 300GB of solid state disk space.  You'll definitely want to increase the CPU as you gain more simultaneous users and RAM depending on your dataset sizes.  Once you create and start the instance:

NOTE: Extra apt install commands are present in the R setup and in the python setup.

```bash
    sudo apt update
    sudo apt upgrade
    sudo apt install build-essential \
        curl \
        git \
        rsync \
        vim \
        unzip \
        wget \
        ca-certificates \
        fonts-roboto \
        fontconfig
    sudo fc-cache -f -v    # cache the fonts
```

Reboot if there are kernel updates (or just to be safe if you don't know.)

### If VM has Hyperdisks

Doing a git checkout with the Hyperdisk mounted may overwrite the symlinked "datasets" and "uploads" directories.

To fix:

`ln -s  /mnt/disks/datastore/datasets /var/www/datasets`
`ln -s  /mnt/disks/datastore/uploads /var/www/uploads`

### MYSQL

`sudo apt install mysql-server`

Follow instructions in [the MySQL setup](./mysql.md) document

### R

Not necessary if you want projectR to run on a Google Cloud Run service (configurable in gear.ini)

Please consult [the R setup](./r_rpy2.md) for packages to install in order to install R and requisite R packages

### RabbitMQ

Not necessary if you want projectR to run in the Apache environment or do not want to setup the RabbitMQ messaging service (configurable in gear.ini)

`sudo mkdir -p /var/log/gEAR_queue`
`sudo chown -R www-data:www-data /var/log/gEAR_queue`

Follow instructions in [the RabbitMQ setup](./rabbitmq.md) document

### Python

Follow instructions in [the Python setup](./python.md) document

### APACHE

`sudo apt install apache2 apache2-dev`

Follow instructions in [the Apache setup](./apache.md) document

### Plotly image export

Install Chromium for kaleido so "Download image" and display thumbnails work. Follow [Plotly image export](./plotly_image_export.md).

### Sass

Used for changing theme colors between portal flavors (e.g. `www/css/gear-theme-purple.scss`).

Ruby Sass (`ruby-sass`) is end-of-life and may not be packaged on current Ubuntu releases. Install Dart Sass 1.79.0 or newer instead, either via npm or a standalone release from <https://github.com/sass/dart-sass/releases>:

```bash
# Remove ruby sass if present, to avoid conflicts with dart-sass
sudo apt remove ruby-sass

# Change VER as needed
VER=1.83.4
cd /tmp
wget https://github.com/sass/dart-sass/releases/download/$VER/dart-sass-$VER-linux-x64.tar.gz
sudo tar -xzf dart-sass-$VER-linux-x64.tar.gz -C /opt
sudo ln -sf /opt/dart-sass/sass /usr/local/bin/sass
hash -r

sass www/css/gear-theme-purple.scss www/css/gear-theme-purple.css
```

Select the desired portal's brand colors in the SCSS before compiling; if multiple
sets are uncommented, the last set takes precedence. The compiled theme also exports
the primary color and a contrast-safe text color used for selected results.

### gEAR portal

```bash
cd ~/git
git clone https://github.com/IGS/gEAR.git
cd /var
sudo rm -rf www && sudo ln -s ~/git/gEAR/www
```

### Systemd Services

More information about these services can be found in [systemd.md](./systemd.md), [RabbitMQ Consumers](../services/rabbitmq_consumers.md) and [Spatial Panel](../services/spatial.md).

```bash
cd ~/git/gEAR/systemd
sudo cp *.target *.service *.slice /etc/systemd/system/

# Replace <gear_root> with the gEAR root on this server.
# spatial-panel.service also has <domain url> and a Python path to check by hand.
sudo sed -i 's|<gear_root>|/path/to/gEAR|g' /etc/systemd/system/*-consumer@.service /etc/systemd/system/spatial-panel.service

# Reload systemd
sudo systemctl daemon-reload

# Enable and start every consumer at once. gear-consumers.target pulls in each
# consumer's target, which pulls in its workers.
sudo systemctl enable gear-consumers.target
sudo systemctl start gear-consumers.target

sudo systemctl enable spatial-panel.service
sudo systemctl start spatial-panel.service

# <Add more services here>
```

To deploy code or unit file changes to the consumers later, follow [Rolling Out Updates](../services/rabbitmq_consumers.md#rolling-out-updates).

### Cron Jobs

`bin/generate_intro_stats.py` pre-computes the statistics, "What's new" lists and study list shown on the intro page (`intro.html`), the featured studies on the home page (`index.html`) and the Studies page (`studies.html`). It reads the database using the credentials in `gear.ini` and writes `www/stats_archive/intro_stats.json`. Until it has run, those sections stay empty or hidden, so run it once by hand and then schedule it daily.

Add the job to the crontab of the account that owns the gEAR checkout, since it needs write access to `www/stats_archive/`. If gEAR uses a virtual environment, put that environment's `python3` in front of the script instead of relying on the `#!/usr/bin/env python3` line, because cron runs with a minimal `PATH`.

```bash
# Cron is installed and running by default on Ubuntu; check with:
systemctl status cron

# Run it once now and confirm it writes the file
cd ~/git/gEAR/bin && ./generate_intro_stats.py

# Add a daily 3:15 AM run to your crontab without opening an editor
(crontab -l 2>/dev/null; echo "15 3 * * * cd $HOME/git/gEAR/bin && ./generate_intro_stats.py >> $HOME/generate_intro_stats.log 2>&1") | crontab -

# Confirm the entry
crontab -l
```

To edit or remove the entry later, use `crontab -e`. The log file collects the script's output and any warnings, such as a `www/landing/<name>/study.json` it skipped because it was invalid.

Rerun the script by hand after adding or editing a study's `study.json`, rather than waiting for the nightly run.

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

### E-mail support

The portal uses e-mail only when users create a new account or use the "forgot password" link
functionality. Here are the steps for setting this up (assuming running on a GCP node):

- $ sudo apt install sendmail
- Go to your Google Account Settings
- Make sure that 2-step verification is turned ON
- Go to your Google Account App Passwords page (https://myaccount.google.com/apppasswords)
- Give your app a custom name, like 'gear-emailer' and click Create
- Copy the 16-character password that appears in the model
- Put your email and this password in your portal's gear.ini file


### Undocumented currently


- What to do for helpdesk config?
