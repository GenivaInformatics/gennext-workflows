# Somatic WGS CNV for the mixed-run scheduling test

`WGS_hg38` selects hg38 canonical chromosomes (1–22, X, Y, M). It is a
somatic DNA kit; the germline sample in the same manifest still uses `wgs`
and does not participate in somatic post-processing.

## Computation

The reference builder recognizes whole-chromosome BED intervals against the exact
FASTA index. Explicit `--mode wgs` is also available. WGS intersects that BED with
CNVkit's accessible-genome regions, annotates and splits those intervals into
10 kb average bins, creates no antitarget bins, and uses `--no-edge` consistently
in reference and fix. The bin size is an engineering starting point, configurable
on the builder; changing it requires matching cohort/batch parameters and a new
reference cache. It is not a validated analytical resolution claim.

Both the main cohort-CNV stage and post-processing use the healthy D4 cohort.
WGS post-processing copies its already-built reference only after exact target-BED
comparison; it never pools tumors into their own normal reference. Existing panel
post-processing retains hybrid binning and its previous sample-pooling behavior.
A mismatched cached WGS mode, BED or bin size fails rather than silently reusing
an incompatible reference.

CNVkit guidance: https://cnvkit.readthedocs.io/en/master/pipeline.html and
https://cnvkit.readthedocs.io/en/v0.8.5/nonhybrid.html .

## Integration requirements

- Publish builder `0.1.2-wgs` using `images/cnv-cohort-builder/Dockerfile.wgs`.
  It inherits the pinned dependency image used in validation.
- Apply `per-sample-cnv-cohort-analysis.yaml`, `post-process-cnv.yaml`, and
  `per-sample-cnv-analysis.yaml` together.
- Merge-patch `argo/gennext-references` with `deploy/wgs-hg38-reference-patch.yaml`.
  The current loader splits kit names at `_`, so the map key is `WGS`, not `WGS_hg38`.
- Deploy the companion orchestrator change: a single `WGS_hg38` tumor is sufficient
  for CNV because its normal reference is external. Other kits keep a minimum of two.
- Only after those pieces are deployed, set this kit's `post_processing_steps` to
  `ARRAY['cnv']`. Do not enable HRD or variant-counter as part of this change.
- The full healthy-reference build happens lazily on first use. Its CPU, memory,
  storage and wall-clock cost must be observed before calling the first sample ready.

The orchestrator companion branch starts from the currently deployed
`feature/post-processing` branch. It does not merge or replace PRs #10–#13 or
finish their remaining mixed-process split.

## Validation

Run unit tests:

```sh
python3 -m unittest discover -s images/cnv-cohort-builder -p 'test_*.py' -v
```

Run isolated Argo validation on acurare-01:

```sh
tests/wgs-cnv/create-validation-configmap.sh
argo submit --context default -n argo tests/wgs-cnv/argo-validation.yaml
```

This uses a synthetic three-normal D4 cohort and a tumor with a known twofold
gain. It exercises real D4 extraction, reference build, cache reuse/mismatch
rejection, CNVkit fix, CBS segmentation and calling. It then generates bins from
the node's actual hg38 FASTA and probes its actual healthy-cohort D4. References
are mounted read-only; all results are under the validation container's `/tmp`.
No patient run is created and no backend completion callbacks are sent.

Raw BAM coverage, full-scale healthy-reference construction, downstream clinical
annotation and patient-specific CNV accuracy are separate validation steps. Supply
the selected BAM/VCF or FASTQ inputs to perform sample-level validation.

## Germline progress

The active frontend polls orchestrator sample status every 10 seconds for all
active process types. The production 0.1.5 orchestrator reads process-type-specific
expected pod totals, counts Argo pod completion events, and returns per-sample
progress. Backend callbacks separately maintain sample and aggregate run status.

Three completed WGS FASTQ workflows had 25 pods, matching their annotation. A
completed WGS BAM workflow had 24 successful pods and no failed/retry pod; its
annotation was 23. This change corrects the WGS BAM total to 24. Pod totals are
step-count estimates, not a fraction of elapsed runtime.

Validation recorded 2026-09-25:
- Six focused Python regression tests passed.
- Real CNVkit 0.9.13 smoke test passed locally and in Argo: 98 synthetic bins,
  recovered gain-vs-baseline log2 contrast 0.9600 for an expected contrast of 1.
- Argo workflow `wgs-cnv-validation-qgg44` succeeded on `acurare-01`:
  292,401 real hg38 bins across 24 contigs; 50 healthy-cohort tracks; all 20 probed
  D4 intervals returned data. The kit BED itself includes 25 contigs; accessible
  regions/bin generation yielded 24. No claim of mitochondrial CNV support.
- Argo lint against live dependencies and syntax checks for six shell templates passed.
- The companion orchestrator's full `go test -count=1 ./...` passed.

Production templates, images, kit post-processing settings and reference maps have
not been changed by this validation. Only the isolated validation Workflow and its
code ConfigMap were submitted. The full patient/mixed-run test still needs input paths.
