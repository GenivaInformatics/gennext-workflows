# NanOnCT Panel v1.0 (Nanodigmbio)

Hybrid-capture somatic DNA panel on MGI, prepared with the NadPrep DNA Library
Preparation Kit v2 and BMI (duplex UMI) adapters. About 73 genes, 13 MSI
markers, the TERT promoter and fusion introns for ALK, BRAF, EGFR, MET, NTRK2
and ROS1 (hg38).

## Pieces

| Piece | Where | Value |
| --- | --- | --- |
| Kit row | `enrichment_kits` in gennext-db and gennext-db-dev | `NanOnCT_Panel_v1.0`, `som-test-tumor`, post-processing `{cnv}` |
| SNV regions | kit row `filename`, under `hg38/regions/` | `NanOnCT_Panel_v1.0_target_hg38.bed` (vendor target BED) |
| CNV regions | `argo/gennext-references` via `deploy/nanoct-hg38-reference-patch.yaml` | covered BED, design id `NanOnCT` |
| CNVkit reference | `hg38/other/cnvkit/NanOnCT_Panel_v1.0_covered_hg38.valid/` | built from the WGS healthy cohort by `ensure-kit-reference` on first use |
| UMI | `umi-detect-template`, image `gennext-umi-detect:0.1.3` | scheme `nadprep_bmi`: 64 fixed families, `5M2S+T` on R1 and R2, duplex |

The SNV step pads the target BED itself, so the target file is used as shipped.
CNV uses the covered (probe) footprint.

## Files for `hg38/regions/`

| File | SHA-256 source | Notes |
| --- | --- | --- |
| `NanOnCT_Panel_v1.0_target_hg38.bed` | vendor file, unchanged | 895 intervals, 258,288 bp |
| `NanOnCT_Panel_v1.0_covered_hg38.bed` | vendor file with CRLF removed | 832 intervals, 374,693 bp; the `-source` key |
| `NanOnCT_Panel_v1.0_covered_hg38.valid.bed` | bed-validator 0.2.1, `--canonical-only` | 3 columns, same coordinates; hg38 confirmed |

Copy all three to `/runspace/gennext/data/ref/hg38/regions/` on every node
that will run the kit, before admitting NanOnCT samples there.

## Apply

```bash
kubectl -n argo patch configmap gennext-references --type merge \
  --patch-file deploy/nanoct-hg38-reference-patch.yaml
kubectl apply -f tasks/base_tasks/umi-detect.yaml
```

The reference patch is partial ConfigMap data; never apply it as a whole
ConfigMap. The first NanOnCT CNV run builds the kit's CNVkit reference, which
takes longer than later runs.

## Open

- No MSI BED is configured (`msi_bed_filename` is empty), so MSI is not scored
  for this kit yet, although the panel carries MSI markers.
- fgbio grouping and duplex consensus have not been checked on NadPrep data.
