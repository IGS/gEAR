FROM jupyter/minimal-notebook:python-3.11

USER root

RUN mkdir -p /opt/gear_starters
COPY notebooks/python_notebook_template.ipynb /opt/gear_starters/python_notebook_template.ipynb
RUN chown -R ${NB_UID}:${NB_GID} /opt/gear_starters

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

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

USER ${NB_UID}
