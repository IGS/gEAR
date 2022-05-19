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

    $ sudo apt install hdf5-helpers hdf5-tools libhdf5-dev zlib1g-dev libblas-dev liblapack-dev libxml2-dev cmake apache2 apache2-dev

    $ ./pip3 install --upgrade pip

    $ ./pip3 install \
      anndata==0.7.8 \
      biocode==0.10.0 \
      biopython==1.79 \
      dash-bio==1.0.2 \
      Flask==2.1.0 \
      Flask-RESTful==0.3.9 \
      h5py==3.6.0 \
      itsdangerous==2.1.2 \
      jupyter==1.0.0 \
      kaleido==0.2.1 \
      llvmlite==0.38.0 \
      mod-wsgi==4.9.0 \
      mysql-connector-python==8.0.28 \
      numba==0.55.1 \
      numexpr==2.8.1 \
      numpy==1.21.5 \
      opencv-python==4.5.5.64 \
      pandas==1.4.1 \
      Pillow==9.0.1 \
      plotly==5.6.0 \
      python-dotenv==0.20.0 \
      requests==2.27.1 \
      rpy2==3.5.1 \
      scanpy==1.8.2 \
      scanpy[louvain]==1.8.2 \
      scikit-learn==1.0.2 \
      scipy==1.8.0 \
      SQLAlchemy==1.4.32 \
      xlrd==1.2.0

SAdkins note - Installing dash-bio (0.6.1) on the devel server by itself upgraded plotly to 5.3.1.  So we need to install dash==1.2.1 and plotly==4.14.3 if they are not this.

    $ sudo mkdir /opt/bin
    $ sudo ln -s /opt/Python-${PYTHONV}/bin/python3 /opt/bin/

Scanpy (or dependencies) assumes it can write in several directories which the web server won't be able to write to by default, so this can be fixed with:

    $ cd /opt/Python-${PYTHONV}/lib/python${PYTHON_MINORV}/site-packages/scanpy
    $ find ./ -name __pycache__ -exec chmod 777 {} \;

Diffxpy (v0.7.4) is used in the multigene curator, but their version does not allow for free-ordering of conditions (for volcano plots).  I forked their code (adkinsrs/diffxpy), made the fix, and installed
    $ /opt/Python-${PYTHONV}/bin/python3 -m pip install git+https://github.com/adkinsrs/diffxpy.git@b2ebeb0fb7c6c215d51264cd258edf9d013ff021

The MulticoreTSNE module currently fails with cmake 3.22.0 or greater.  I have a pending pull request to fix this but until then:
    $ /opt/Python-${PYTHONV}/bin/python3 -m pip install git+https://github.com/jorvis/Multicore-TSNE.git@68325753c4ab9758e3d242719cd4845d751a4b6c

