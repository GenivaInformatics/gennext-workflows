    # Migrating from `parabricks-mutect2-tumor-only` to a sharded GATK CPU step

## Why we're moving

Two real-world miss cases driving this change:

- **CDKN2A** — sample 103027-2, chr9:21974697 TAACTATTC>T, 4.1% AF. Manta caught it as a small indel; Parabricks Mutect2 (4.3.1-1 *and* 4.7.0-1) did not, regardless of `--mutect-low-memory`, `--pileup-detection`, or any LOD/active-region knob.
- **STAT3** — sample 103174-5, chr17:42322474 C>T, 4.7% AF. No secondary caller catches it (SNP, not SV), so we cannot force-call it. Parabricks misses it because its `--pileup-detection` enforces a hard-coded 10% AF floor for SNV pileup injection, and the underlying GATK4 flag `--pileup-detection-snp-alt-threshold` is **not** exposed through `--mutectcaller-options`.

Lifting the caller from Parabricks to **GATK 4.6 CPU** makes the full GATK4 sensitivity surface available — `--pileup-detection-snp-alt-threshold`, `--pileup-detection-enable-indel-pileup-calling`, `--use-pdhmm`, plus the Broad-recommended FFPE workflow (`--f1r2-tar-gz` → `LearnReadOrientationModel` → `FilterMutectCalls --orientation-bias-artifact-priors`). Mutect2 itself is single-threaded, so we recover wallclock by sharding the panel and scattering Mutect2 across N parallel JVMs inside the same pod.

## Architecture

Single Argo step, single container, internal scatter-gather. From Argo's perspective the contract is identical to the old `parabricks-mutect2-tumor-only` step (same parameter names, same input BAM mount, single output VCF). The image swaps from `nvcr.io/nvidia/clara/clara-parabricks` to `broadinstitute/gatk` and the `--gpus all` requirement disappears.

Inside the container:

1. `gatk SplitIntervals --scatter-count N --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION` — produces N balanced interval lists from the panel BED. `BALANCING_WITHOUT_INTERVAL_SUBDIVISION` keeps each capture region whole, so reads spanning interval boundaries don't get split across shards.
2. `xargs -P N` over the shard ids — N parallel `gatk Mutect2` JVMs, each scoped to its own interval list, emitting a per-shard VCF *and* a per-shard F1R2 tar.
3. `gatk MergeVcfs` + `gatk MergeMutectStats` + `gatk LearnReadOrientationModel` — consolidate the N partial outputs. All three are required: `FilterMutectCalls` reads the merged stats, and the orientation model needs the F1R2 tars from every shard pooled together to converge.
4. `gatk FilterMutectCalls --orientation-bias-artifact-priors` — applies the FFPE-aware filter pass, plus the standard Mutect2 filters (germline, contamination, strand bias, etc.).

The Broad image (`broadinstitute/gatk:4.6.2.0`) ships GATK, bcftools 1.13, samtools 1.13, tabix, bgzip, and htsfile, so the entire pipeline runs in one image without any add-on tooling.

## Container script (drop into the Argo step)

```bash
#!/usr/bin/env bash
set -euo pipefail

REF=/mnt/ref/${REF_FA}
BAM=/mnt/input/bam/${BAM_BASENAME}
PANEL=/mnt/input/bed/${PANEL_BASENAME}
OUT_VCF=/mnt/output/vcf/${OUT_VCF_BASENAME}
SAMPLE_ID=${SAMPLE_ID}
SHARDS=${SHARDS:-10}
HMM_THREADS=${HMM_THREADS:-2}

WORK=/mnt/scratch/${SAMPLE_ID}
mkdir -p "$WORK"/{shards,vcfs,f1r2,logs}

# 1. SplitIntervals
gatk SplitIntervals \
  -R "$REF" -L "$PANEL" \
  --scatter-count "$SHARDS" \
  --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
  -O "$WORK/shards"

# 2. Parallel Mutect2
run_shard() {
  local shard="$1"
  gatk --java-options "-Xmx4g" Mutect2 \
    -R "$REF" \
    -I "$BAM" \
    -L "$WORK/shards/${shard}-scattered.interval_list" \
    --tumor-sample "$SAMPLE_ID" \
    --pileup-detection \
    --pileup-detection-snp-alt-threshold 0.02 \
    --pileup-detection-enable-indel-pileup-calling \
    --pileup-detection-indel-alt-threshold 0.02 \
    --pileup-detection-absolute-alt-depth 3 \
    --use-pdhmm \
    --f1r2-tar-gz "$WORK/f1r2/${shard}.tar.gz" \
    --native-pair-hmm-threads "$HMM_THREADS" \
    -O "$WORK/vcfs/${shard}.vcf.gz" \
    > "$WORK/logs/${shard}.log" 2>&1
}
export -f run_shard
export REF BAM WORK SAMPLE_ID HMM_THREADS

ls "$WORK/shards"/*-scattered.interval_list \
  | sed 's|.*/||;s|-scattered.interval_list||' \
  | xargs -P "$SHARDS" -I {} bash -c 'run_shard "$@"' _ {}

# 3. Merge
VCF_ARGS=();   for f in "$WORK"/vcfs/*.vcf.gz;       do VCF_ARGS+=(-I  "$f"); done
STATS_ARGS=(); for f in "$WORK"/vcfs/*.vcf.gz.stats; do STATS_ARGS+=(-stats "$f"); done
F1R2_ARGS=();  for f in "$WORK"/f1r2/*.tar.gz;       do F1R2_ARGS+=(-I "$f"); done

gatk MergeVcfs           "${VCF_ARGS[@]}"   -O "$WORK/${SAMPLE_ID}.unfiltered.vcf.gz"
gatk MergeMutectStats    "${STATS_ARGS[@]}" -O "$WORK/${SAMPLE_ID}.unfiltered.vcf.gz.stats"
gatk --java-options "-Xmx8g -XX:+UseParallelGC" LearnReadOrientationModel \
                         "${F1R2_ARGS[@]}"  -O "$WORK/${SAMPLE_ID}.priors.tar.gz"

# 4. Filter (FFPE orientation-bias filter included)
gatk FilterMutectCalls \
  -R "$REF" \
  -V "$WORK/${SAMPLE_ID}.unfiltered.vcf.gz" \
  --stats "$WORK/${SAMPLE_ID}.unfiltered.vcf.gz.stats" \
  --orientation-bias-artifact-priors "$WORK/${SAMPLE_ID}.priors.tar.gz" \
  -O "$OUT_VCF"
```

Two implementation notes worth carrying forward:

- **`-Xmx8g -XX:+UseParallelGC` on `LearnReadOrientationModel` is mandatory.** With default G1GC and the default heap, the EM step crashed with a JVM SIGSEGV on the higher-coverage 103174-5 sample. ParallelGC plus an 8 GiB heap is stable; this is the only tool in the chain that needs a non-default JVM config.
- **`--native-pair-hmm-threads 2` per shard, not more.** The active-region traversal is single-threaded and is the bottleneck; pair-HMM threads beyond 2 produce diminishing returns and over-subscribe the scheduler when N shards run concurrently. Per-shard CPU sits at ~110% (1 active worker + ~1 HMM worker), so a 10-shard pod uses ~13 cores in steady state.

## Argo WorkflowTemplate (drop-in over `parabricks-mutect2-tumor-only`)

```yaml
inputs:
  parameters:
    - name: input-bam
    - name: sample-id
    - name: genome-ref
    - name: regions-file
    - name: output-vcf
    - { name: shards,      value: "10" }
    - { name: hmm-threads, value: "2" }
volumes:
  - { name: ref,    hostPath: { path: "{{=sprig.osDir(inputs.parameters['genome-ref'])}}",   type: Directory } }
  - { name: bam,    hostPath: { path: "{{=sprig.osDir(inputs.parameters['input-bam'])}}",    type: Directory } }
  - { name: bed,    hostPath: { path: "{{=sprig.osDir(inputs.parameters['regions-file'])}}", type: Directory } }
  - { name: vcf,    hostPath: { path: "{{=sprig.osDir(inputs.parameters['output-vcf'])}}",   type: DirectoryOrCreate } }
  - { name: scratch, emptyDir: { sizeLimit: "10Gi" } }
script:
  image: broadinstitute/gatk:4.6.2.0
  resources:
    requests: { cpu: "16",  memory: "24Gi" }
    limits:   { cpu: "20",  memory: "32Gi" }
  env:
    - { name: REF_FA,           value: "{{=sprig.osBase(inputs.parameters['genome-ref'])}}" }
    - { name: BAM_BASENAME,     value: "{{=sprig.osBase(inputs.parameters['input-bam'])}}" }
    - { name: PANEL_BASENAME,   value: "{{=sprig.osBase(inputs.parameters['regions-file'])}}" }
    - { name: OUT_VCF_BASENAME, value: "{{=sprig.osBase(inputs.parameters['output-vcf'])}}" }
    - { name: SAMPLE_ID,        value: "{{inputs.parameters.sample-id}}" }
    - { name: SHARDS,           value: "{{inputs.parameters.shards}}" }
    - { name: HMM_THREADS,      value: "{{inputs.parameters.hmm-threads}}" }
  command: ["bash", "-c"]
  source: |
    # paste the container script from above here verbatim
volumeMounts:
  - { name: ref,     mountPath: /mnt/ref,        readOnly: true }
  - { name: bam,     mountPath: /mnt/input/bam,  readOnly: true }
  - { name: bed,     mountPath: /mnt/input/bed,  readOnly: true }
  - { name: vcf,     mountPath: /mnt/output/vcf }
  - { name: scratch, mountPath: /mnt/scratch }
```

The differences from the existing Parabricks template are limited to:

- **Image**: `broadinstitute/gatk:4.6.2.0` instead of `nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1`.
- **No GPU**: drop `nvidia.com/gpu` from `resources.limits` and remove any GPU-related toleration / node selector. Schedule on CPU pool.
- **`--in-tumor-recal-file` is not needed**: the BAM was already recal-applied via `pbrun applybqsr` upstream. Mutect2 reads the recalibrated BAM directly.
- **Scratch volume**: an `emptyDir` of ~10 GiB for shard intermediates. Sized to fit ~50 MB per shard × 10 shards plus interval lists and logs.
- **Inputs/outputs identical**: same parameters, same paths — the rest of the DAG (annotation, downstream filtering, MAF generation) doesn't need to change.

## Resource budget and runtime

Validated on the 64-core / 128 GiB host:

| | wallclock | CPU cores used (steady state) | RSS |
|---|---|---|---|
| Parabricks 4.3.1-1 (current production) | ~100 s | 1 GPU (A4000) + ~3 CPU | ~6 GiB |
| Single-thread GATK CPU `--pileup-detection` | ~22 min | ~2 cores | ~2 GiB |
| Sharded GATK CPU, 10 shards × 2 HMM threads (this design) | **~3-5 min** | ~13-15 cores | ~10-12 GiB |

Pod sizing: request 16 CPU / 24 GiB, limit 20 CPU / 32 GiB. The CPU request gives the scheduler enough room for the 10 active-region workers plus auxiliary HMM threads; the memory request covers ten 4 GiB-heap shards plus the LearnReadOrientationModel 8 GiB peak and ample buffer.

If wallclock matters more than node-packing, raise `SHARDS` to 16 or 20 — the panel has 12,542 intervals and SplitIntervals balances cleanly down to ~600 intervals/shard. Beyond ~20 shards the per-JVM startup and merge-stage costs start to dominate and you'll see diminishing returns.

## FFPE rationale (why these specific flags)

Three classes of sensitivity / specificity adjustment, each driven by a real failure mode:

- `--pileup-detection-snp-alt-threshold 0.02` and `--pileup-detection-enable-indel-pileup-calling --pileup-detection-indel-alt-threshold 0.02 --use-pdhmm` — bypass the de Bruijn assembler for low-AF SNVs and indels. Without these, the CDKN2A and STAT3 cases are dropped before they reach the haplotype-likelihood stage.
- `--pileup-detection-absolute-alt-depth 3` — guards against over-emission at very low depth. The relative threshold of 2% is meaningless when there are only two reads at the locus; this requires at least three alt reads in absolute terms.
- `--f1r2-tar-gz` + `LearnReadOrientationModel` + `FilterMutectCalls --orientation-bias-artifact-priors` — Broad's official FFPE filter. Cytosine deamination during formalin fixation produces sample-specific C>T / G>A artifacts at low AF that look identical to real low-AF somatic events on AF alone. The orientation model learns the per-sample F1R2/F2R1 imbalance and tags artifact-pattern reads. On 103174-5 the STAT3 C>T receives ROQ=19 (orientation-bias quality, lower = more suspect) but still PASSes — the filter is doing its job: it scrutinized a real C>T while keeping it.

## Migration path / validation checklist

- [ ] Run the new template on samples 103027-2 and 103174-5; confirm both target variants are PASS in the filtered output. Expected TLODs: 29.54 (CDKN2A) and 27.41 (STAT3).
- [ ] Run on 5–10 historical samples; diff the filtered VCFs against the production Parabricks output. Flag any PASS calls present in Parabricks but absent in the new pipeline (regression risk), and any new PASS calls (the rescue we're paying for).
- [ ] Build a panel of normals from 30–50 archived FFPE-normals if available, replacing the gnomAD-only default — captures recurrent FFPE artefact loci that the orientation filter alone misses.
- [ ] Once confidence is high, retire the `parabricks-mutect2-tumor-only` template.

## Outputs

The new step emits:

- `${output-vcf}` — filtered VCF (PASS / non-PASS), drop-in for whatever currently consumes the Parabricks output.
- `${output-vcf}.unfiltered.vcf.gz` (optional, can keep in scratch) — pre-filter for downstream debugging.
- `${output-vcf}.stats` — Mutect2 call stats, used by FilterMutectCalls and useful for review.
- `${sample-id}.priors.tar.gz` (optional) — orientation-bias priors, useful for QC and re-filtering.

Decide whether to expose the auxiliary files as workflow outputs or keep them only in the pod scratch; the filtered VCF is the only one downstream stages need.
