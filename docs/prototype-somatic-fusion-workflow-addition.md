# RNA-seq fusion workflow: architecture decisions

## Context

Designing an RNA-seq fusion detection workflow for FFPE tumor samples. This issue captures the architectural decisions and the reasoning behind them, so reviewers don't re-litigate the same questions and so future maintainers understand what's intentional vs. accidental.

## Workflow topology

```
Input FASTQs
    ↓
fastp (adapter trim; UMI extraction when applicable)
    ↓
SortMeRNA (rRNA quantification + depletion)
    ↓
    ├──→ STAR (Arriba params) ──────→ Arriba ────────┐
    │                                                 ├──→ fusion-report ──→ FusionInspector ──→ MultiQC
    └──→ STAR (STAR-Fusion params) ──→ STAR-Fusion ──┘                                            ↑
                            │                                                                      │
                            └──→ Picard MarkDup ──→ RustQC ───────────────────────────────────────┘
                                 (QC sidecar)        (Qualimap, RSeQC + TIN, Preseq,
                                                     dupRadar, featureCounts, Samtools)
```

## Decisions

### Base pipeline: nf-core/rnafusion 4.x

Provides the wiring (fastp, two-STAR fork, Arriba, STAR-Fusion, fusion-report, FusionInspector, Picard MarkDup, MultiQC) out of the box. Reference build is ~24h on first run. FusionCatcher is skipped to avoid the COSMIC account requirement and reference-build burden.

We bolt on SortMeRNA upstream of the STAR fork and RustQC downstream of MarkDup. UMI handling is deferred to a separate branch (see UMI section).

### Two STAR alignments, not one

Arriba and STAR-Fusion need different chimeric parameter sets. Running STAR once with shared params is technically possible (OICR-GSI ships a clinical workflow doing this with `--chimOutType 'WithinBAM HardClip Junctions'`) but costs Arriba sensitivity at the margins, which is the wrong tradeoff for FFPE.

| Tool | `--chimOutType` | `--chimSegmentMin` | Notes |
|------|-----------------|---------------------|-------|
| Arriba | `WithinBAM HardClip` | 10 | Liberal scoring, recovers weak chimerics |
| STAR-Fusion | `Junctions WithinBAM SoftClip` | 12 | Tighter scoring, `peOverlapNbasesMin 12` |

nf-core/rnafusion maintainers chose the two-STAR split deliberately (issue #436). We follow that lead. The unified-STAR approach is acceptable for development/benchmarking runs where iteration speed matters more than the last 5% of recall, but production stays on two STARs.

Compute cost: ~2x alignment time vs. unified, but alignment is not the bottleneck for clinical-throughput volumes. Worth it for the sensitivity floor on FFPE.

### SortMeRNA: rRNA quantification + depletion

FFPE libraries fail rRNA depletion at variable rates — bad lots of Twist or KAPA RiboErase produce libraries that are 30-60% rRNA, and there's no way to tell from FastQC alone. SortMeRNA classifies reads against curated rRNA databases (SILVA / Rfam) and provides:

- **`% reads classified as rRNA`** — the QC metric that goes to MultiQC.
- **Depleted FASTQs** — rRNA-removed reads passed downstream to STAR.

Configured to both quantify and filter. Rationale for filtering:

- rRNA reads cannot contribute supporting evidence to any fusion call, so removal is non-destructive for the fusion-calling goal.
- Recovers effective coverage budget for the on-target fraction (a "50M-read library that's 40% rRNA" becomes a 30M-read library on-target, and STAR's resource needs scale with the smaller input).
- Original raw FASTQs are preserved upstream of fastp/SortMeRNA, so the operation is reversible.

Sits between fastp and the STAR fork so both STAR runs ingest the same depleted FASTQs.

### FusionInspector for IGV visualization

FusionInspector is **downstream of the callers**, not parallel to them. Its input is a fusion candidate list (`GeneA--GeneB` pairs, one per line) plus the trimmed FASTQs; for each candidate it builds a fusion contig (gene A + gene B joined at the breakpoint), realigns reads to those contigs, and produces a consolidated BAM ready for IGV.

Two invocation modes exist:

- **Integrated** with STAR-Fusion via `--FusionInspector inspect|validate` — auto-pulls STAR-Fusion's own calls only.
- **Standalone** with `--fusions <list>` — accepts an external candidate list.

We use the standalone mode, fed by **fusion-report's merged candidate list**. This is what nf-core/rnafusion does and is the right choice here because it pulls candidates from *both* Arriba and STAR-Fusion (and FusionCatcher if/when re-enabled), not just STAR-Fusion. Anything either caller flags gets visualized; the integrated mode would silently drop Arriba-only calls.

Outputs we care about:

- `finspector.consolidated.cSorted.bam` + `.bai` — the IGV target.
- `finspector.fa` + `.fai` — fusion contig reference for IGV.
- `finspector.bed` / `finspector.gtf` — exon track for fused contigs.
- `finspector.fusion_inspector_web.html` — interactive web view.
- `finspector.igv.xml` — pre-built IGV session file.

For FFPE, override the defaults that drop low-expression candidates:

```
--min_FFPM 0
--min_junction_reads 1
--vis
--examine_coding_effect
```

Whole-gene contig approach is preferred over exact-base custom contigs because FusionInspector's index includes decoy genes and other candidates, so ambiguous reads are pulled away appropriately. Custom contigs lose this and overcall on paralogous regions.

### NVIDIA Parabricks: rejected for this stage

`nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1` exposes only a stripped `pbrun starfusion` that takes `Chimeric.out.junction` as input — no FASTQ entry point, no `--FusionInspector` flag. FusionInspector itself is CPU-only.

Parabricks pays off for high-throughput `pbrun rna_fq2bam` GPU alignment, but since FusionInspector runs on CPU per sample regardless, the throughput ceiling is set elsewhere. Use `trinityctat/starfusion:latest` for the FusionInspector step.

Revisit if/when sample volume makes GPU alignment worthwhile and FusionInspector gets sharded across nodes.

### Picard MarkDuplicates: QC-only sidecar

MarkDup runs on the STAR-Fusion-tuned BAM and feeds RustQC + MultiQC. It does **not** feed back into Arriba or STAR-Fusion.

- Arriba does its own coordinate-based dedup internally (`duplicates` filter).
- STAR-Fusion calls primarily off `Chimeric.out.junction`, not the BAM. The maintainer's position (STAR-Fusion issue #134 and related threads) is that STAR-Fusion's evidence model is robust to duplicates without requiring pre-emptive removal.
- Coordinate-based dedup on RNA-seq is unreliable without UMIs — highly expressed transcripts produce many genuine fragments with identical 5′ coords, and Picard can't tell those from PCR duplicates.

This matches nf-core/rnafusion's convention.

### RustQC: consolidated post-alignment QC

Runs on the STAR-Fusion BAM after MarkDup. Single-pass tool covering what would otherwise be 5+ separate steps:

| RustQC module | Replaces | Notes |
|---------------|----------|-------|
| dupRadar | R dupRadar package | PCR-vs-natural duplication diagnostic; emits MultiQC custom-content files |
| featureCounts | Subread featureCounts | Gene-level quantification |
| RSeQC suite | RSeQC tools | Includes `infer_experiment`, `read_distribution`, `junction_saturation`, `inner_distance`, `bam_stat`, `read_duplication`, **TIN** |
| Qualimap RNA-seq | Qualimap | 3'/5' bias, gene-body coverage |
| Preseq | Preseq `lc_extrap` | Library complexity / saturation curve |
| Samtools stats | samtools | flagstat / idxstats / stats |

TIN (Transcript Integrity Number) is particularly relevant for FFPE — alignment-based RNA degradation metric, complementary to RIN. Vanilla `tin.py` is impractically slow (~10h on a 186M-read BAM); RustQC computes it as part of the single-pass run.

All RustQC outputs are MultiQC-compatible by design.

### fastp: trim only, no FASTQ-stage dedup

Configuration for the non-UMI path:

```
fastp \
  --in1 R1.fq.gz --in2 R2.fq.gz \
  --out1 R1.trim.fq.gz --out2 R2.trim.fq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 25 \
  --json sample.fastp.json --html sample.fastp.html
```

`--dedup` is **not** enabled. It performs sequence-hash deduplication at the FASTQ stage, which (a) is not UMI-aware even when `--umi` is set (UMIs live in read names, not the sequence fastp hashes), and (b) compounds the same RNA-seq coordinate-dedup problem more aggressively than Picard MarkDup, since it lacks even mate-pair coordinate context.

### UMI handling: deferred branch

Current workflow assumes non-UMI library prep. If FFPE protocols move to UMI'd panels, the dedup picture changes:

1. **fastp UMI extraction** — `--umi --umi_loc=read1 --umi_len=N --umi_skip=N` per kit spec. Replaces the separate `umi_tools extract` step; UMI ends up appended to the read name.
2. STAR alignment as usual (both Arriba- and STAR-Fusion-tuned passes).
3. `umi_tools dedup` on each BAM — sets `BAM_FDUP` correctly using UMI + position.
4. Arriba runs with the `-u` flag to honor external dup marks (Arriba is otherwise oblivious to UMIs and dedups solely on coordinates).
5. Picard MarkDup becomes redundant; `umi_tools dedup` produces the dup-rate metrics.

This is a separate workflow branch, not a rewiring of the existing path.

## Single MultiQC report

All steps are wired to emit MultiQC-compatible outputs into a single report:

- fastp — native module
- SortMeRNA — native module
- STAR (×2) — native module
- Picard MarkDup + CollectRnaSeqMetrics — native module
- RustQC (Qualimap, RSeQC, Preseq, dupRadar, featureCounts, Samtools, TIN) — native modules + custom content for dupRadar
- Arriba / STAR-Fusion / FusionInspector / fusion-report — surfaced via custom content

## Open questions

- [ ] UMI prep adoption timeline — drives the priority of the umi_tools branch.
- [ ] Transcript-level quantification (Salmon / RSEM) — RustQC's featureCounts covers gene-level. Add transcript-level if isoform-aware analysis or fusion-partner-isoform readouts are needed.
- [ ] Sample identity check (Somalier) — clinical safety net for cohort-scale processing. Cheap to add but not yet decided.
- [ ] Whether to publish a single-STAR benchmarking run on our FFPE cohort to quantify the sensitivity delta vs. two-STAR. Useful for periodic sanity-checking the tradeoff.
- [ ] Argo `WorkflowTemplate` design for the FusionInspector step — fits as a downstream template after fusion calling.

## References

- nf-core/rnafusion: https://nf-co.re/rnafusion
- nf-core/rnafusion two-STAR rationale: https://github.com/nf-core/rnafusion/issues/436
- OICR-GSI unified-STAR workflow: https://github.com/oicr-gsi/arriba
- Arriba algorithm docs: https://arriba.readthedocs.io/en/latest/internal-algorithm/
- Arriba UMI handling: https://arriba.readthedocs.io/en/latest/current-limitations/
- FusionInspector: https://github.com/FusionInspector/FusionInspector/wiki
- STAR-Fusion duplicate-removal discussion: https://github.com/STAR-Fusion/STAR-Fusion/issues/134
- SortMeRNA: https://github.com/biocore/sortmerna
- RustQC: https://seqeralabs.github.io/RustQC/
- nf-core/rnaseq (reference QC bundle): https://nf-co.re/rnaseq
- MultiQC supported tools: https://docs.seqera.io/multiqc/modules
- Trinity CTAT image: `trinityctat/starfusion:latest`