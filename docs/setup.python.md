## Overview

These are instructions for setting up the Python environment for Python
environment, tested on Ubuntu 18.04 and 20.04.  Why not virtualenv or any
other isolated python environment?  You could, but for a committed
webserver that's an unnecessary layer.  Also, maybe I'm old-school, but
fixed paths have worked fine for decades.

    $ sudo apt install libffi-dev libsqlite-dev libsqlite3-dev libhdf5-dev
    $ sudo apt install libreadline-gplv2-dev libncursesw5-dev libssl-dev python3-tk tk-dev libgdbm-dev libc6-dev libbz2-dev lzma liblzma-dev

    $ sudo mkdir /opt/Python-3.7.3
    $ sudo chown $USER /opt/Python-3.7.3
    $ cd /tmp
    $ wget https://www.python.org/ftp/python/3.7.3/Python-3.7.3.tar.xz
    $ tar -xf Python-3.7.3.tar.xz
    $ cd Python-3.7.3/
    $ ./configure --prefix=/opt/Python-3.7.3 --enable-optimizations --enable-shared
    $ make
    $ make install

    $ cd /opt/Python-3.7.3/bin
    $ sudo ln -s /opt/Python-3.7.3/lib/libpython3.7m.so.1.0 /usr/lib/

    $ sudo apt install hdf5-helpers hdf5-tools libhdf5-dev zlib1g-dev libblas-dev liblapack-dev libxml2-dev cmake apache2 apache2-dev

    $ ./pip3 install --upgrade pip
    $ ./pip3 install h5py==2.10.0 scanpy==1.5.1 anndata==0.7.3 pandas==1.0.4 numba==0.50.0 xlrd==1.2.0 numpy==1.18.5 scipy==1.4.1 scikit-learn==0.23.1 jupyter mysql-connector-python==8.0.20 requests plotly==4.14.3 kaleido==0.2.1 louvain MulticoreTSNE Pillow biopython==1.76 biocode==0.9.0 python-dotenv==0.14.0 Flask==1.1.2 SQLAlchemy==1.2.12 Flask-RESTful==0.3.8 mod-wsgi==4.6.5 opencv-python==4.1.1.26 pathlib==1.0.1 dash-bio==0.6.1

    $ sudo mkdir /opt/bin
    $ sudo ln -s /opt/Python-3.7.3/bin/python3 /opt/bin/

Scanpy (or dependencies) assumes it can write in several directories which the web server won't be able to write to by default, so this can be fixed with:

    $ cd /opt/Python-3.7.3/lib/python3.7/site-packages/scanpy
    $ find ./ -name __pycache__ -exec chmod 777 {} \;
