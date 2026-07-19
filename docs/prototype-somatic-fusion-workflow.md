# RNA-seq fusion workflow: architecture decisions

## Context

Designing an RNA-seq fusion detection workflow for FFPE tumor samples. This issue captures the architectural decisions and the reasoning behind them, so reviewers don't re-litigate the same questions and so future maintainers understand what's intentional vs. accidental.

## Workflow topology

```
Input FASTQs
    ↓
fastp (adapter trim, optional UMI handling)
    ↓
    ├──→ STAR (Arriba params) ──→ Arriba ──┐
    │                                       ├──→ fusion-report ──→ FusionInspector ──→ MultiQC
    └──→ STAR (STAR-Fusion params) ──→ STAR-Fusion ──┘                                  ↑
                            │                                                            │
                            └──→ Picard MarkDup ──→ RustQC ─────────────────────────────┘
                                 (QC sidecar)
```

## Decisions

### Base pipeline: nf-core/rnafusion 4.x

Provides the wiring (fastp, two-STAR fork, Arriba, STAR-Fusion, fusion-report, FusionInspector, Picard MarkDup, CollectRnaSeqMetrics, MultiQC) out of the box. Reference build is ~24h on first run. FusionCatcher is skipped to avoid the COSMIC account requirement and reference-build burden.

Gaps we bolt on:
- **RustQC** as a post-pipeline step on the STAR-Fusion-tuned BAM (replaces / supplements the default QC modules with our internal tooling).
- **UMI handling** — not in the pipeline yet; deferred to a separate branch (see UMI section below).

### Two STAR alignments, not one

Arriba and STAR-Fusion need different chimeric parameter sets. Running STAR once with shared params is technically possible (OICR-GSI ships a clinical workflow doing exactly this with `--chimOutType 'WithinBAM HardClip Junctions'`) but costs Arriba sensitivity at the margins, which is the wrong tradeoff for FFPE.

| Tool | `--chimOutType` | `--chimSegmentMin` | Notes |
|------|-----------------|---------------------|-------|
| Arriba | `WithinBAM HardClip` | 10 | Liberal scoring, recovers weak chimerics |
| STAR-Fusion | `Junctions WithinBAM SoftClip` | 12 | Tighter scoring, `peOverlapNbasesMin 12` |

nf-core/rnafusion maintainers chose the two-STAR split deliberately (issue #436). We follow that lead. The unified-STAR approach is acceptable for development/benchmarking runs where iteration speed matters more than the last 5% of recall, but production stays on two STARs.

Compute cost: ~2x alignment time vs. unified, but alignment is not the bottleneck for clinical-throughput volumes. Worth it for the sensitivity floor on FFPE.

### FusionInspector for IGV visualization

FusionInspector builds fusion contigs (gene A + gene B joined at the breakpoint), realigns reads to those contigs, and produces a consolidated BAM ready for IGV. Outputs we care about:

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

MarkDup runs on the STAR-Fusion-tuned BAM and feeds CollectRnaSeqMetrics + MultiQC. It does **not** feed back into Arriba or STAR-Fusion.

- Arriba does its own coordinate-based dedup internally (`duplicates` filter).
- STAR-Fusion calls primarily off `Chimeric.out.junction`, not the BAM.
- Coordinate-based dedup on RNA-seq is unreliable without UMIs — highly expressed transcripts produce many genuine fragments with identical 5′ coords, and Picard can't tell those from PCR dupes.

This matches nf-core/rnafusion's convention.

### UMI handling: deferred branch

Current workflow assumes non-UMI library prep. If FFPE protocols move to UMI'd panels, the dedup picture changes:

1. `umi_tools extract` between fastp and STAR.
2. STAR alignment (both passes).
3. `umi_tools dedup` on each BAM, setting `BAM_FDUP` properly using UMI + position.
4. Arriba runs with `-u` to honor external dup marks (Arriba is otherwise oblivious to UMIs).
5. Picard MarkDup becomes redundant.

This is a separate workflow branch, not a rewiring of the existing path.

### RustQC integration

Runs on the STAR-Fusion BAM after MarkDup. Outputs converge with FusionInspector results into MultiQC. Replaces/supplements default QC modules with our internal tooling — keeps QC outputs aligned with the rest of our reporting.

## Open questions

- [ ] UMI prep adoption timeline — drives the priority of the umi_tools branch.
- [ ] Whether to publish a single-STAR benchmarking run on our FFPE cohort to quantify the sensitivity delta vs. two-STAR. Useful for periodic sanity-checking the tradeoff.
- [ ] Argo `WorkflowTemplate` design for the FusionInspector step — fits as a downstream template after fusion calling.

## References

- nf-core/rnafusion: https://nf-co.re/rnafusion
- nf-core/rnafusion two-STAR rationale: https://github.com/nf-core/rnafusion/issues/436
- OICR-GSI unified-STAR workflow: https://github.com/oicr-gsi/arriba
- Arriba algorithm docs: https://arriba.readthedocs.io/en/latest/internal-algorithm/
- Arriba UMI handling: https://arriba.readthedocs.io/en/latest/current-limitations/
- FusionInspector: https://github.com/FusionInspector/FusionInspector/wiki
- Trinity CTAT image: `trinityctat/starfusion:latest`