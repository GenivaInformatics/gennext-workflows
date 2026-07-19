FROM mambaorg/micromamba:1.5.10

USER root

RUN micromamba install -y -n base -c bioconda -c conda-forge \
        cnvkit=0.9.13 \
        d4tools=0.3.11 \
        numpy \
    && micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH

COPY cnv_cohort_builder.py /usr/local/bin/cnv-cohort-builder
RUN chmod +x /usr/local/bin/cnv-cohort-builder

ENTRYPOINT ["cnv-cohort-builder"]
