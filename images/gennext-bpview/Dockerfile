FROM mambaorg/micromamba:1.5.10

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
        samtools=1.21 \
        pysam=0.22 \
        python=3.12 \
    && micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH

COPY bpview.py /opt/bpview/bpview.py
COPY tests /opt/bpview/tests
RUN python /opt/bpview/tests/test_bpview.py \
    && ln -s /opt/bpview/bpview.py /usr/local/bin/bpview \
    && chmod +x /opt/bpview/bpview.py

# No ENTRYPOINT: the Argo template sets `command:`.
