FROM jupyter/minimal-notebook:python-3.11

USER root

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER ${NB_UID}

RUN pip install --no-cache-dir \
    scanpy \
    anndata \
    h5py \
    numpy \
    scipy \
    pandas \
    matplotlib \
    seaborn \
    scikit-learn \
    notebook \
    jupyterlab

RUN mkdir -p /home/nbuser/gear_starters

COPY notebooks/python_starter.ipynb /home/nbuser/gear_starters/python_starter.ipynb

USER root
RUN chown -R ${NB_UID}:${NB_GID} /home/nbuser/gear_starters
USER ${NB_UID}
