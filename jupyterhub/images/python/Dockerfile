FROM jupyter/minimal-notebook:python-3.11

USER root

RUN mkdir -p /home/nbuser/gear_starters

COPY notebooks/python_notebook_template.ipynb /home/jovyan/gear_starters/python_notebook_template.ipynb

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

USER root
RUN chown -R ${NB_UID}:${NB_GID} /home/jovyan/gear_starters
USER ${NB_UID}
