# Python setup instructions

## Overview

These are instructions for setting up the Python environment, tested on Ubuntu 22.04.  Why not virtualenv or any
other isolated python environment?  You could, but for a committed
webserver that's an unnecessary layer.  Also, maybe I'm old-school, but
fixed paths have worked fine for decades.

```bash
    sudo apt-get -qq update
    sudo DEBIAN_FRONTEND="noninteractive" apt -qq install -y --no-install-recommends \
        apache2-dev \
        libffi-dev \
        libsqlite3-dev \
        libreadline-dev \
        libncursesw5-dev \
        libssl-dev \
        tk-dev \
        libgdbm-dev \
        libc6-dev \
        liblzma-dev \
        libbz2-dev \
        zlib1g-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libblosc-dev \
        libhdf5-dev \
        hdf5-helpers \
        hdf5-tools \
        libblas-dev \
        liblapack-dev \
        libcairo2-dev \
        libxml2-dev \
        libwebp-dev \
        libpcre2-dev \
        libicu-dev \
        libdeflate-dev \
        libssl3 \
        pkg-config \
        llvm \
        apache2 \
        php \
        libapache2-mod-php \
        php-gd \
    sudo apt -qq clean autoclean \
    sudo apt -qq autoremove -y \
    sudo rm -rf /var/lib/apt/lists/*

    export PYTHONV=3.14.4
    export PYTHON_MINORV=3.14

    sudo mkdir /opt/Python-${PYTHONV}
    sudo chown $USER /opt/Python-${PYTHONV}
    cd /tmp
    wget https://www.python.org/ftp/python/${PYTHONV}/Python-${PYTHONV}.tar.xz
    tar -xf Python-${PYTHONV}.tar.xz
    cd Python-${PYTHONV}/
    # LDFLAGS, CPPFLAGS, and PKG_CONFIG_PATH are needed to ensure the build process can find the lzma library and headers,
    # which are required for Python's lzma module.  Without these flags, the build process may fail to find lzma downstream.
    ./configure --prefix=/opt/Python-${PYTHONV} --enable-optimizations --enable-shared \
        LDFLAGS="-L/usr/lib/x86_64-linux-gnu" \
        CPPFLAGS="-I/usr/include" \
        PKG_CONFIG_PATH="/usr/lib/x86_64-linux-gnu/pkgconfig" \
    make install
    cd /tmp
    rm -rf Python-${PYTHONV}*

    export LD_LIBRARY_PATH="/opt/Python-${PYTHON_FULL_VERSION}/lib:${LD_LIBRARY_PATH}"

    cd /opt/Python-${PYTHONV}/bin
    sudo ln -s /opt/Python-${PYTHONV}/lib/libpython${PYTHON_MINORV}.so.1.0 /usr/lib/

    sudo mkdir /opt/bin
    sudo ln -s /opt/Python-${PYTHONV}/bin/python3 /opt/bin/

    sudo apt install r-base r-base-dev hdf5-helpers hdf5-tools libhdf5-dev zlib1g-dev libblas-dev liblapack-dev libxml2-dev cmake apache2 apache2-dev

    # Disable user-installs to guarantee packages go into /opt/Python...
    export PIP_USER=0
```

## pip install option A (reqs file)

Check the requirement.txt file in <git_repo_root>/docker for the latest packages that have been tested locally. They work in a Dockerized ubuntu environment so should be able to work on the VMs. You can run `./pip3 install -r requirements.txt` as a shortcut.

`./pip3 install -r <git_repo_root/docker/requirements.txt`
`./pip3 uninstall dask-expr -y`
`./pip3 install -e <git_repo_root>/lib/`

If the requirements.txt will not install due to a stack depth issue, you can use the `requirements.full.txt` instead, which was made using `pip freeze > requirements.full.txt`. This file contains all versioned scripts so pip does not have to compute the best version for non-mentioned packages.

## pip install option B (manual)


NOTE: Some of the packages will indirectly install dask-expr, which is currently broken for spatialdata, with no intention of fixing. So it is necessary to uninstall dask-expr to avoid issues with spatialdata (<https://github.com/scverse/spatialdata/pull/570>)

NOTE 2: Really try to keep the requirements.txt in sync with the files below.  While you can use the file to do the installation, @jorvis prefers to have them explicitly listed here.

```bash
    ./pip3 install --upgrade pip

    ./pip3 install \
    aiohttp==3.13.5 \
    aiohttp_retry==2.9.1 \
    anndata==0.12.11 \
    biocode==0.10.0 \
    cairosvg==2.7.1 \
    colorcet==3.1.0 \
    dash-bio==1.0.2
    datashader==0.19.0 \
    Flask==3.1.3 \
    Flask-RESTful==0.3.9 \
    gosling==0.3.0 \
    hic2cool==0.8.3 \
    jupyterlab==4.0.5 \
    jupyter==1.0.0 \
    kaleido==0.2.1 \
    leidenalg==0.10.2 \
    legacy-cgi==2.6.4 \   # To handle legacy CGI scripts (cgi was removed in Python 3.13)
    llvmlite==0.47.0 \
    matplotlib==3.10.7 \
    mod-wsgi==5.0.2 \
    more_itertools==11.0.2 \
    mysql-connector-python==8.0.28 \
    numba==0.65.0 \
    numpy==2.4.0 \
    opencv-python==4.5.5.64 \
    openpyxl==3.1.5 \
    pandas==2.3.3 \   # Anndata does not support pandas 3.x.x yet
    panel==1.8.10 \
    Pillow==12.2.0 \
    pika==1.3.2 \
    pims==0.7.0 \ # Need to force the latest version, due to numpy being pretty recent. Relates to spatialdata-io dependency chain
    plotly==6.6.0 \
    pybigwig==0.3.25 \
    python-dotenv==0.20.0 \
    requests==2.31.0 \
    rpy2==3.6.7 \
    scanpy==1.12.1 \
    scipy==1.17.1 \
    seaborn==0.13.2 \
    setuptools<82 \   # need some pkg_resources methods
    spatialdata==0.7.2 \
    spatialdata_io==0.6.0 \
    shadows==0.1a2 \
    tables==3.11.1 \
    watchfiles==1.1.1
    ./pip3 install git+https://github.com/adkinsrs/diffxpy.git@ffd828c280882ca98adc6e42c934625fab0011f6
    ./pip3 uninstall dask-expr -y

```

### Note about editable pip installs

The previous pip installation methods also includes an extra line to install the gEAR "lib" area as an editable install. Updates to the modules in this directory will be hot-loaded without a re-install.  The "setup.py" script is designed to find the "lib" directory itself as a package. This will allow you to run `import gear` without having to append "lib" into the PYTHONPATH.  However, you will still need to append "lib" to the PYTHONPATH if you want to `import geardb` or something on the same level as "setup.py"

## Gotchas

Scanpy (or dependencies like numba) assumes it can write in several directories which the web server won't be able to write to by default, so this can be fixed with:

```bash
    cd /opt/Python-${PYTHONV}/lib/python${PYTHON_MINORV}/site-packages/
    find ./ -name __pycache__ -exec chmod 777 {} \;
```

## Note about Flask

The Flask server is set to run with debug mode as False by default, but by setting the DEBUG environment variable to "True" you can change this (since the value is read from os.environ)
