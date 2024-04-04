FROM continuumio/miniconda3:main as base

WORKDIR /app
RUN ls /home
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda install --yes segmm -c sihanliu
