# TST500C CNV reference

The TST500C pooled CNVkit reference on `gennext-02` is keyed by
`TST500C_manifest_hg38.valid.bed`. The validated BED has SHA-256
`bb64b36e5d1793e64b53a6be90869de956358cdcf11f3d97237425e522cc1689`.
Its `build.json` records `refFlat_mane.txt`, stored in
`hg38/gene_annots/`. The global `hg38.regions.ref-flat` points to a
different file and must not be used to validate this reference.

Install the validated BED at
`/runspace/gennext/data/ref/hg38/regions/TST500C_manifest_hg38.valid.bed`
on each node that will run TST500C CNV. Verify its hash against `build.json`
before admitting work. Merge-patch `argo/gennext-references` with
`deploy/tst500c-hg38-reference-patch.yaml`, then apply
`workflows/per-sample-cnv-cohort-analysis.yaml`. The reference patch is
partial ConfigMap data; applying it as a complete ConfigMap would remove
unrelated reference keys.
