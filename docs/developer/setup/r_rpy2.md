# Setting up R, projectR, and the rpy2 module

Rather than using the OS' package manager to install R, we'll use the steps/script below to do it so we can more carefully control versioning.

## Prerequesites to install via apt

```bash
sudo apt -qq update
sudo DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
  gfortran \
  libreadline-dev \
  xorg-dev \
  libbz2-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libfontconfig1-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libglpk-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libwebp-dev \
  libgit2-dev \
  libuv1-dev \
  pandoc \
  tzdata
sudo apt -qq clean autoclean
sudo apt -qq autoremove -y
sudo rm -rf /var/lib/apt/lists/*
```

## Installing R

Run this to install R, Bioconductor, and projectR:

from the gEAR root:

```bash
cd ./docker/
# Install R and packages
sudo sh ./install_R.sh

export R_HOME="/usr/local/lib/R"
export LD_LIBRARY_PATH="/usr/local/lib/R/lib:${LD_LIBRARY_PATH}"

```

### If install fails

I have observed "pak" installation failure as it relies on an "apt update" internally and will fail if that has errors.

Some possible fixes (based on errors encountered)

```bash
# Import the missing Google Cloud GPG key so APT can verify the Google repository
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo gpg --dearmor -o /usr/share/keyrings/google-cloud-keyring.gpg
echo "deb [signed-by=/usr/share/keyrings/google-cloud-keyring.gpg] https://packages.cloud.google.com/apt google-cloud-ops-agent-jammy-all main" | sudo tee /etc/apt/sources.list.d/google-cloud-ops-agent.list

# Remove or Disable Broken RabbitMQ Repositories
sudo rm -f /etc/apt/sources.list.d/*rabbitmq*.list
```

Another fix is to ensure `options(pkg.sysreqs = FALSE)` in install_packages.R.  It is true for Docker installs but should be false for server installs.

Run the install after this.

### Fix the libR.conf file

To ensure R's shared libraries are found create a file "libR.conf" in /etc/ld.so.conf.d and add the following contents:

```bash
# libR default configuation
/usr/local/lib/R/lib
```

Then run `sudo ldconfig` to cache the shared libraries.  You can confirm shared libaries with `sudo ldconfig -v`

## Installing rpy2

This will be done in the Python installation doc steps.
