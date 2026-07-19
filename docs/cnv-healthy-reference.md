# CNV Healthy Reference Set — WGS D4 Storage and Per-Kit Reference Build

This document describes how the germline "healthy" cohort is stored and how a per-kit CNVkit pooled reference is generated from it for each clinical cancer enrichment kit (CGP panels, etc.).

The pipeline targets **somatic CNV calling on FFPE tumor samples** against a stable, kit-agnostic source of normals.

## 1. Why this layout

CNVkit's `pooled_reference.cnn` is bin-locked to one specific target+antitarget BED. Reusing it across kits is unsupported (see CNVkit issues #85, #525, #864 and the same constraint in GATK CNV / PureCN / ichorCNA). The community SOP — which `nf-core/createpanelrefs` formalizes — is "one source cohort, regenerate the PoN per kit."

We adopt that pattern, with two refinements:

1. **WGS germlines as the universal source.** Every clinical kit's targets fall inside a WGS cohort, so one cohort serves all current and future panels. WES germlines would only cover kits whose targets sit inside the exome capture.
2. **D4 depth tracks instead of BAMs.** For a CNV-only health set we don't need alignments — only per-base depth. D4 (mosdepth `--d4`) compresses ~20–50× vs. BAM, supports fast random-access aggregation over arbitrary BEDs, and is what the per-kit build job actually consumes.

## 2. Storage layout

```
${REF_BASE_DIR}/healthy-cohort-wgs/
├── manifest.tsv                # one row per sample (schema below)
├── d4/
│   ├── HC-001.d4
│   ├── HC-002.d4
│   └── ...
├── kits/
│   ├── agilent-cgp-dna-rna/
│   │   ├── targets.bed         # cnvkit.py target output
│   │   ├── antitargets.bed     # cnvkit.py antitarget output
│   │   ├── pooled_reference.cnn
│   │   ├── coverage/           # per-sample synthetic .cnn (kept for audit)
│   │   │   ├── HC-001.targetcoverage.cnn
│   │   │   ├── HC-001.antitargetcoverage.cnn
│   │   │   └── ...
│   │   └── build.json          # build metadata (date, cohort-hash, tool versions)
│   ├── twist-comprehensive-exome/
│   └── ...
└── flat/
    └── (per-kit no-cohort fallback references)
```

**`d4/`** is the source of truth. All kit-specific artifacts under `kits/` are reproducible from `d4/` + the kit BED + the reference genome. They can be regenerated and should be regenerated when the cohort changes.

## 3. Manifest schema

`manifest.tsv` (tab-separated, header included):

| Column | Description |
|---|---|
| `sample-id` | Internal cohort ID (e.g. `HC-001`) |
| `sex` | `M`, `F`, `unknown` — used for chrX/chrY handling and per-sex sub-references |
| `source-tissue` | `blood`, `saliva`, `buffy-coat` — FFPE-derived germlines flagged separately (see §6) |
| `platform` | Sequencer model (e.g. `NovaSeq 6000`, `NovaSeq X Plus`) |
| `chemistry` | Library prep + flowcell (e.g. `TruSeq DNA PCR-Free, S4 v1.5`) |
| `read-length` | e.g. `2x150` |
| `mean-coverage` | Mean WGS coverage from `mosdepth` summary |
| `date-sequenced` | YYYY-MM-DD |
| `consent-id` | Consent / IRB tracking ID |
| `notes` | Free text (ancestry if known, exclusions, etc.) |

The build job reads this manifest to (a) pick which samples participate in a given build, (b) split by sex when constructing sex-specific references, and (c) record cohort composition in `build.json`.

## 4. Adding a sample to the cohort

Inputs: a germline WGS BAM aligned to the same reference genome we use for tumor analysis (currently hg38).

```bash
# One-time per sample. Run on a node with the BAM accessible.
mosdepth \
    --threads 4 \
    --no-per-base \
    --fast-mode \
    --d4 \
    ${REF_BASE_DIR}/healthy-cohort-wgs/d4/${SAMPLE_ID} \
    ${BAM}
```

This produces `${SAMPLE_ID}.d4` plus `*.mosdepth.summary.txt` (used to populate `mean-coverage` in the manifest).

Then append a row to `manifest.tsv`. Drop the BAM once D4 is verified, per the cohort retention policy.

**Inclusion criteria.** Aim for ≥30 samples (CNVkit recommends ≥10 minimum, ≥30 for stable noise estimates). Match the production tumor population on:

- Sequencer model and chemistry — different platforms have different GC and edge biases.
- Read length.
- Ancestry diversity proportional to the patient population, where consent permits.
- Sex balance (or build sex-specific sub-references; see §5).

## 5. Building a per-kit reference

This is an **offline batch job**, not part of the per-run post-process workflow. Run it once when onboarding a new kit, and re-run when the cohort composition changes meaningfully or when the kit BED is revised.

### Inputs

- `${KIT_BED}` — kit's target regions BED (vendor or custom)
- `${REF_GENOME}` — genome FASTA (must match what the d4 files were generated against)
- `${REF_FLAT}` — refFlat for gene annotation
- The cohort `d4/` directory and `manifest.tsv`

### Steps

**Step 1 — Build the bin grid for this kit.** Standard CNVkit `target` / `antitarget` / `access`, identical to what `cnvkit-prepare-references` does in `workflows/post-process-cnv.yaml`:

```bash
cnvkit.py access ${REF_GENOME} -o access.bed

cnvkit.py target ${KIT_BED} \
    --annotate ${REF_FLAT} \
    --split \
    -o targets.bed

cnvkit.py antitarget targets.bed \
    -g access.bed \
    -o antitargets.bed
```

**Step 2 — Synthesize per-sample `.cnn` files from D4.** For each cohort sample, aggregate D4 depth over `targets.bed` and `antitargets.bed` and emit CNVkit-compatible coverage files.

`.targetcoverage.cnn` / `.antitargetcoverage.cnn` are tab-separated with the columns:

```
chromosome  start  end  gene  log2  depth
```

For each bin, `depth` is the mean per-base depth in `[start, end)`, and `log2` is `log2(depth / median_bin_depth)` (CNVkit's convention; the median is computed within target and antitarget pools separately). `gene` is taken from the annotated `targets.bed`; antitarget rows use `Antitarget`.

This is a small custom step — CNVkit has no native D4 importer. Recommended implementation: a Python script using `pyd4` (or `d4tools view --bed`) to extract mean depth per interval, emitted as TSV.

**Step 3 — Pool into a kit reference.** Standard CNVkit:

```bash
cnvkit.py reference \
    coverage/*.targetcoverage.cnn coverage/*.antitargetcoverage.cnn \
    -f ${REF_GENOME} \
    -o pooled_reference.cnn
```

Add `--sample-sex` and split by sex from the manifest if building sex-specific references. CNVkit also accepts `-y` to assume male samples have one X.

**Step 4 — Record build metadata.** Write `build.json`:

```json
{
  "kit": "agilent-cgp-dna-rna",
  "kit-bed-sha256": "...",
  "kit-bed-path": "/ref/regions/agilent_cgp.bed",
  "ref-genome": "hg38",
  "cohort-samples": ["HC-001", "HC-002", "..."],
  "cohort-hash": "sha256 of sorted sample-ids",
  "n-samples": 42,
  "sex-split": {"M": 21, "F": 21},
  "cnvkit-version": "0.9.13",
  "built-at": "2026-04-27T10:00:00Z",
  "built-by": "build-cnv-reference-job/abc123"
}
```

This is what `batch-load-and-prepare` should validate against at run time — if the kit's `build.json.kit-bed-sha256` doesn't match the BED currently configured for the kit, fail loudly rather than running with a stale reference.

## 6. FFPE caveat

The cohort above is **blood/saliva-derived germline WGS**. Tumor samples in production are **FFPE**. FFPE introduces:

- More fragmentation and shorter usable insert sizes.
- Amplified GC bias (AT-rich regions deplete disproportionately).
- Higher and less uniform coverage noise floor.

A reference built only from non-FFPE normals will under-model the FFPE noise floor — `pooled_reference.cnn`'s `spread` column will be tighter than the tumor data warrants, producing **more false-positive low-amplitude CNV calls** in the FFPE tumor.

Mitigation options, in order of preference:

1. **Maintain a parallel FFPE-normal sub-cohort** (FFPE-derived germline material — adjacent normal tissue, FFPE buffy-coat, etc.) when available. Build `kits/{kit}/pooled_reference_ffpe.cnn` from it and route FFPE tumor samples to that reference. This is the cleanest fix but depends on availability of FFPE-normal material.
2. **Augment the WGS cohort with FFPE pseudo-normals** — FFPE tumor samples with regions of known cancer drivers masked, treated as normals for noise modeling only. Suitable when matched FFPE-normals don't exist. Requires careful curation and a documented mask BED.
3. **Tune downstream filtering**, not the reference. Tighten the `.cnr` log2 thresholds and segment-confidence cutoffs in `cnvkit-post-process` for FFPE samples. Cheapest, but doesn't fix the underlying noise mismatch.
4. **Flat reference fallback** (`cnvkit.py reference -f genome.fa -t targets.bed -a antitargets.bed`) for kits with neither blood nor FFPE normals available. Stored under `flat/`. Loses the cohort noise model but is portable; acceptable as a sanity-check baseline, not as primary clinical output.

The current recommendation is to start with the blood/saliva WGS cohort, deploy it for clinical use with downstream filter tuning (option 3), and begin collecting FFPE-derived normals in parallel toward option 1.

## 7. Integration with the existing workflow

`workflows/post-process-cnv.yaml` currently builds the pooled reference inline at runtime from the current batch's tumor samples (lines 100–109). With this design:

- `cnvkit-prepare-references` (lines 164–223) and `cnvkit-build-reference` (lines 278–306) move out of the per-run DAG.
- `batch-load-and-prepare` is extended to look up `${REF_BASE_DIR}/healthy-cohort-wgs/kits/{enrichment-kit}/` and emit `pooled-reference-cnn`, `targets-bed`, `antitargets-bed` as outputs, validated against the kit BED's sha256.
- `cnvkit-coverage` (lines 230–273) reads `targets.bed` / `antitargets.bed` from the kit reference directory rather than `/scratch/references/`.
- `cnvkit-analyze` (`per-sample-cnv-analysis.yaml` lines 237–241) reads `pooled_reference.cnn` from the kit reference directory.

Net effect: per-run CNV becomes coverage + fix + segment + call + downstream, with the reference treated as a static input keyed on `enrichment-kit`. Smaller batches stop producing misleading results because the reference is no longer derived from the batch.

## 8. Maintenance

- **Cohort refresh.** When ≥5 new samples are added or platform/chemistry shifts, re-run the per-kit build job for all active kits. Old `pooled_reference.cnn` files are kept under `kits/{kit}/archive/{date}/` for audit.
- **Kit BED revision.** Vendors occasionally release updated BEDs. A new BED ⇒ a full rebuild for that kit. The `kit-bed-sha256` in `build.json` is the trigger for a stale-reference failure at run time.
- **Genome reference change.** All `d4/` files are tied to the alignment reference. A genome change requires regenerating D4 from the source alignments (kept in cold storage per consent policy) or from re-aligned BAMs.
- **Sample retirement.** If a cohort sample is later found to carry a pathogenic CNV, remove it from the manifest, archive its `.d4`, and rebuild affected kits.

## References

- CNVkit pipeline docs — https://cnvkit.readthedocs.io/en/stable/pipeline.html
- nf-core/createpanelrefs — https://github.com/nf-core/createpanelrefs
- mosdepth (D4 output) — https://github.com/brentp/mosdepth
- D4 format paper — https://www.biorxiv.org/content/10.1101/2020.10.23.352567v1.full
- pyd4 — https://github.com/38/d4-format/tree/main/pyd4
