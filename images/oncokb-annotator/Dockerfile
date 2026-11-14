# syntax=docker/dockerfile:1.7
#
# Upstream oncokb-annotator, unmodified, plus our wrapper scripts.
#
# OncoKB publishes no container image, so we build upstream's source at a
# pinned commit. Nothing under /app/oncokb is patched; the only additions are
# the three wrapper scripts copied in at the end. See README.md.
#
#   0.2.0-upstream-v3.4.2  upstream tag v3.4.2 (1cbdfb2), same code and
#                          Python packages as 0.1.7-unfiltered
#   0.3.0-upstream-v4.1.0  upstream tag v4.1.0 (77521cd); the MAF wrapper
#                          expects v4's mixed-mode MAF header
FROM python:3.12.12-alpine3.22@sha256:848ba4413eb897e225159b8fc1b02094576cbae4aa73fc13142608ae2c8c0e32

ARG ONCOKB_ANNOTATOR_VERSION=v4.1.0
ARG ONCOKB_ANNOTATOR_COMMIT=77521cd2fc568422d87d84b73c24ca65b670a8e4

LABEL org.opencontainers.image.source="https://github.com/GenivaInformatics/gennext-workflows" \
      org.opencontainers.image.description="oncokb-annotator ${ONCOKB_ANNOTATOR_VERSION} (unmodified) with gennext wrapper scripts" \
      io.oncokb.annotator.version="${ONCOKB_ANNOTATOR_VERSION}" \
      io.oncokb.annotator.commit="${ONCOKB_ANNOTATOR_COMMIT}"

ADD https://github.com/oncokb/oncokb-annotator.git#${ONCOKB_ANNOTATOR_COMMIT} /app/oncokb
WORKDIR /app/oncokb

# Upstream's own requirement files, with the unpinned transitive packages held
# at the versions 0.1.7-unfiltered shipped. kiwisolver has no wheel for this
# Python, so g++ is needed only while it builds.
COPY constraints.txt /tmp/constraints.txt
RUN apk add --no-cache --virtual .build-deps g++ \
    && pip install --no-cache-dir -c /tmp/constraints.txt \
        -r requirements/common.txt -r requirements/pip3.txt \
    && apk del .build-deps \
    && rm /tmp/constraints.txt

COPY filtered_maf_annotator.py filtered_sv_annotator.py backfill_filtered_maf.py /app/oncokb/
