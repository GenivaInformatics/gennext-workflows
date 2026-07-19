FROM mambaorg/micromamba:1.5.10

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
        samtools=1.21 \
        pysam=0.22 \
        python=3.12 \
    && micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH

COPY extract_reads.py /usr/local/bin/extract-reads-for-fusion-inspector
COPY propagate_dup_flags.py /usr/local/bin/propagate-fi-dup-flags
RUN chmod +x /usr/local/bin/extract-reads-for-fusion-inspector \
              /usr/local/bin/propagate-fi-dup-flags

# No ENTRYPOINT: pods set `command:` explicitly to choose between the two
# tools (extract-reads-for-fusion-inspector or propagate-fi-dup-flags).
