FROM quay.io/jupyter/r-notebook:latest

USER root

RUN mkdir -p /opt/gear_starters
COPY notebooks/r_notebook_template.ipynb /opt/gear_starters/r_notebook_template.ipynb
RUN chown -R ${NB_UID}:${NB_GID} /opt/gear_starters

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        git \
        curl \
        libhdf5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libpng-dev \
        zlib1g-dev \
        libcairo2-dev \
        libjpeg-dev \
        libtiff-dev \
        libfreetype6-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libglpk-dev \
        libgmp-dev \
        libarpack2-dev \
        libblas-dev \
        liblapack-dev \
        python3-dev \
        python3-venv \
        python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER ${NB_UID}

RUN mamba install -y \
    -c conda-forge \
    -c bioconda \
    r-seurat \
    r-cairo \
    cairo \
    pango \
    fontconfig \
    freetype \
    harfbuzz \
    fribidi \
    bioconductor-zellkonverter \
    bioconductor-singlecellexperiment \
    bioconductor-scater \
    bioconductor-scran \
    && mamba clean -afy

RUN R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org/')"
#RUN R -e "BiocManager::install(c('zellkonverter', 'SingleCellExperiment', 'scater', 'scran'), ask=FALSE, update=FALSE)"
RUN R -e "install.packages(c('igraph', 'png', 'reticulate'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages(c('Seurat'), repos='https://cloud.r-project.org/')"
