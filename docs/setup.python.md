# Python setup instructions

## Overview

These are instructions for setting up the Python environment for Python
environment, tested on Ubuntu 22.04.  Why not virtualenv or any
other isolated python environment?  You could, but for a committed
webserver that's an unnecessary layer.  Also, maybe I'm old-school, but
fixed paths have worked fine for decades.

    $ sudo apt install libffi-dev libsqlite-dev libsqlite3-dev libhdf5-dev
    $ sudo apt install libncursesw5-dev libssl-dev python3-tk tk-dev libgdbm-dev libc6-dev libbz2-dev lzma liblzma-dev

    $ export PYTHONV=3.10.4
    $ export PYTHON_MINORV=3.10
    $ sudo mkdir /opt/Python-${PYTHONV}
    $ sudo chown $USER /opt/Python-${PYTHONV}
    $ cd /tmp
    $ wget https://www.python.org/ftp/python/${PYTHONV}/Python-${PYTHONV}.tar.xz
    $ tar -xf Python-${PYTHONV}.tar.xz
    $ cd Python-${PYTHONV}/
    $ ./configure --prefix=/opt/Python-${PYTHONV} --enable-optimizations --enable-shared
    $ make
    $ make install

    $ cd /opt/Python-${PYTHONV}/bin
    $ sudo ln -s /opt/Python-${PYTHONV}/lib/libpython${PYTHON_MINORV}.so.1.0 /usr/lib/

    $ sudo apt install r-base r-base-dev hdf5-helpers hdf5-tools libhdf5-dev zlib1g-dev libblas-dev liblapack-dev libxml2-dev cmake apache2 apache2-dev

## pip install option A (reqs file)

Check the requirement.txt file in <git_repo_root>/docker for the latest packages that have been tested locally. They work in a Dockerized ubuntu environment so should be able to work on the VMs. You can run `./pip3 install -r requirements.txt` as a shortcut.

``./pip3 install -r <git_repo_root/docker/requirements.txt`

## pip install option B (manual)

    $ ./pip3 install --upgrade pip

    $ ./pip3 install \
      aiohttp==3.8.3 \
      anndata==0.10.6 \
      biocode==0.10.0 \
      biopython==1.79 \
      cairosvg==2.7.1 \
      dash-bio==1.0.2 \
      Flask==3.0.0 \
      Flask-RESTful==0.3.9 \
      gunicorn \
      h5py==3.10.0 \
      jupyterlab==4.0.5 \
      jupyter==1.0.0 \
      kaleido==0.2.1 \
      leidenalg==0.10.2 \
      llvmlite==0.41.1 \
      matplotlib==3.9.0 \
      mod-wsgi==4.9.4 \
      more_itertools==9.0.0 \
      mysql-connector-python==8.0.20 \
      numba==0.58.1 \
      numexpr==2.8.4 \
      numpy==1.26.4 \
      opencv-python==4.5.5.64 \
      openpyxl==3.1.5 \
      pandas==2.2.1 \
      Pillow==10.2.0 \
      pika==1.3.1 \
      plotly==5.6.0 \
      python-dotenv==0.20.0 \
      requests==2.31.0 \
      rpy2==3.5.16 \
      sanic \
      scanpy==1.10.1 \
      scikit-learn==1.0.2 \
      scipy==1.11.04 \
      seaborn==0.13.2 \
      shadows==0.1a0 \
      tables==3.9.2 \
      xlrd==1.2.0
    $ sudo mkdir /opt/bin
    $ sudo ln -s /opt/Python-${PYTHONV}/bin/python3 /opt/bin/

# Gotchas

Scanpy (or dependencies like numba) assumes it can write in several directories which the web server won't be able to write to by default, so this can be fixed with:

    $ cd /opt/Python-${PYTHONV}/lib/python${PYTHON_MINORV}/site-packages/scanpy
    $ find ./ -name __pycache__ -exec chmod 777 {} \;

NOTE: Installing custom version of diffxpy that is based on the latest commit on the main branch (at the time). It does not have a release tag, but fixes a NumPy bug occurs with older diffxpy commits and newer numpy releases.

    $ /opt/Python-${PYTHONV}/bin/python3 -m pip install git+https://github.com/theislab/diffxpy.git@7609ea935936e3739fc4c71b75c8ee8ca57f51ea

The MulticoreTSNE module currently fails with cmake 3.22.0 or greater.  I have a pending pull request to fix this but until then:

    $ /opt/Python-${PYTHONV}/bin/python3 -m pip install git+https://github.com/jorvis/Multicore-TSNE.git@68325753c4ab9758e3d242719cd4845d751a4b6c

## Note about Flask

The Flask server is set to run with debug mode as False by default, but by setting the DEBUG environment variable to "True" you can change this (since the value is read from os.environ)
