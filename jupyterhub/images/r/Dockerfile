FROM jupyter/r-notebook:latest

USER root

RUN mkdir -p /opt/gear_starters
COPY notebooks/r_notebook_template.ipynb /opt/gear_starters/r_notebook_template.ipynb
RUN chown -R ${NB_UID}:${NB_GID} /opt/gear_starters

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libhdf5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER ${NB_UID}

RUN R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('zellkonverter', 'SingleCellExperiment', 'scater', 'scran'), ask=FALSE, update=FALSE)"
RUN R -e "install.packages(c('Seurat'), repos='https://cloud.r-project.org/')"
