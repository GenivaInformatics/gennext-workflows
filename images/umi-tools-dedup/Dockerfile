FROM mambaorg/micromamba:1.5.10

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
        umi_tools=1.1.6 \
        samtools=1.21 \
        python=3.12 \
    && micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH

# No ENTRYPOINT: pods set `command:` explicitly. Both umi_tools and
# samtools are on PATH.
