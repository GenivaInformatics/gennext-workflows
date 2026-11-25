# Arriba's protein-domain annotation (Pfam), for fusion-events' kinase_* columns. Taken from the
# Arriba image the pipeline runs, so both read the same file.
FROM quay.io/biocontainers/arriba:2.5.1--h87b9561_0@sha256:702a3e8db8a9f6e227c0154c460f54fc2dc4c46f16d9e3f53b82edf43f05e881 AS arriba

FROM mambaorg/micromamba:1.5.10

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
        samtools=1.21 \
        pysam=0.22 \
        python=3.12 \
    && micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH

COPY --from=arriba /usr/local/var/lib/arriba/protein_domains_hg38_GRCh38_v2.5.1.gff3 /opt/bpview/ref/
# OncoKB curated fusion list, a pinned snapshot (alteration names): the "exact pair" test of fusion-events.
COPY ref/oncokb_curated_fusions.tsv /opt/bpview/ref/
COPY bpview.py fusion_events.py /opt/bpview/
COPY tests /opt/bpview/tests
RUN python /opt/bpview/tests/test_bpview.py \
    && python /opt/bpview/tests/test_fusion_events.py \
    && ln -s /opt/bpview/bpview.py /usr/local/bin/bpview \
    && ln -s /opt/bpview/fusion_events.py /usr/local/bin/fusion-events \
    && chmod +x /opt/bpview/bpview.py /opt/bpview/fusion_events.py

# No ENTRYPOINT: the Argo template sets `command:`.
