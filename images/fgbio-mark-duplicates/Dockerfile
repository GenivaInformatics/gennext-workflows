FROM mambaorg/micromamba:1.5.10

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
        fgbio=2.5.21 \
        samtools=1.21 \
        pysam=0.22 \
        python=3.12 \
        openjdk=21 \
    && micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH

COPY move_umi_to_rx.py /usr/local/bin/move-umi-to-rx
RUN chmod +x /usr/local/bin/move-umi-to-rx
