# Live WorkflowTemplates vs `main` (2026-11-14)

Refs #33. This compares the `argo` namespace today (read-only `kubectl get`) with `main` at 1e1fcd6, the `main` of the time of the comparison. It updates `snapshots/live-2026-11-13/DRIFT.md`; since that snapshot, the rollout applies brought `duckdb-template` and `manta-template` back in line with `main`.

**`main` has moved on since:** at 170b3bf, #45 changes `tasks/base_tasks/manta.yaml` only. Against 170b3bf the same export gives 127 matches and 19 drifts; the extra one is `manta-template`, repo-ahead until #45 is applied live. `oncokb-o2c-dev` is a temporary dev template created after the comparison and is not counted.

**Method**
- Templates are compared by `spec`, as in #28. A repo template is any `WorkflowTemplate` document outside `snapshots/`.
- **Use** counts workflows still on the cluster (14 Sep to 14 Nov 2026, 568 workflows) that reach a template, directly or through other templates' `templateRef`s. A count of 0 means unused for two months, not unused forever.
- **Direction:** the live spec is checked against every past version of its repo file.
  - It matches an old version: the repo moved on and live was never re-applied ("repo ahead").
  - It matches no version: live carries changes that were never committed ("live ahead").
  - The diffs in `drift/` run repo → live.

| Status | Count |
|---|---|
| Same spec as `main` | 128 |
| Drift: name in the repo, spec differs | 18 |
| Live only | 51 |
| Repo only (not on the cluster) | 15 |

This directory changes no existing template and applies nothing.

## Drift (18)
Runs = workflows in the last two months that reached the template.

| Template | Runs | Direction | What differs (repo → live) | Proposed fix |
|---|---|---|---|---|
| `copy-temp-to-output` | 533 | live ahead | Cleanup uses `find … -exec rm -rf` instead of `rm -r /mnt/temp/*` (handles empty dirs and dotfiles) | Commit live |
| `fastqc-template` | 533 | live ahead | Live adds workflow-level `arguments` (a node-name default), `entrypoint`, `nodeSelector` and `/opt/gennext` hostPath `volumes`. The `fastqc` template itself is identical; the extra fields only apply when the template is submitted on its own | Restore repo: drop the workflow-level fields |
| `parabricks-bqsr-template` | 482 | live ahead | Maps dbSNP155 to dbSNP157 for older submissions; fails early if the dbSNP file or its index is missing | Commit live |
| `filter-msisensor-pro-template` | 417 | repo ahead | The repo adds a second template, `filter-msisensor-pro-top10pct`; live is the 2025-09-10 version | Apply repo (gate:live-template); nothing calls the new template yet |
| `extract-reads-for-fusion-inspector-template` | 147 | live ahead | Image `extract-fusion-reads:0.3.3` (repo 0.3.0); new input `artifact-filter-mode` (default `filter`) | Commit live |
| `somatic-rna-fastq-v2` | 147 | live ahead | Passes `artifact-filter-mode: passthrough` to `extract-fusion-reads` | Commit live |
| `somatic-rna-bam-v2` | 147 | live ahead | Arriba and bpview GTF from `ref-paths` (`$.genes.annotation`) instead of a hard-coded MANE path | Commit live |
| `copy-inputs-to-scratch` | 62 | live ahead | Script runs with `set -eu` | Commit live |
| `exomedepth-template` | 62 | live ahead | Trailing newline at the end of the R script only | Commit live (cosmetic) |
| `post-process-cnv` | 9 | live ahead | New `prepare-cnv-regions` step: validates each sample's CNV BED with `gennext-bed-validator:0.2.1` and hands `batch-load-and-prepare` a per-run config. The token comes from the `gennext-cluster-secrets` Secret, not a literal | Commit live |
| `somatic-fusion-fastq-gpu-v2` | 0 | live ahead | Same `artifact-filter-mode: passthrough` as `somatic-rna-fastq-v2` | Commit live |
| `somatic-fusion-fastq-v2` | 0 | diverged | Live has no `extract-fusion-reads` step: FusionInspector reads fastp's FASTQs and depends on `fusion-report`. The repo has the extract step. Live was last updated 2026-05-13, the repo 2026-07-19 | Decide with the FusionInspector retirement (wave 2); until then, apply repo |
| `msi-tmb-coverage-workflow-template` | 0 | live ahead | `collect-tmb-data` gets `filter-multihit: true` and `multihit-method: poisson`; `cosmic-vcf` is no longer passed | Commit live, after confirming with the TMB owner |
| `wrapper-workflow-prs` | 0 | live ahead | Passes `genotype-file-name: <product-id>.tsv.gz` | Commit live (PRS owner to confirm) |
| `genotype-conditional-workflow-template` | 0 | live ahead | Trailing newline only | Commit live (cosmetic) |
| `annotsv-to-oncokb-template` | 0 | repo ahead | The repo adds `annotsv-to-oncokb-on-node`; live is the 2024-04-04 version | Apply repo, or retire both if unused |
| `trim-galore-template` | 0 | repo ahead | The repo caps `--cores` at 16 and quotes `--basename`; live is the 2023-10-09 version | Apply repo |
| `bed-split-workflow` | 0 | diverged | Kebab-case parameters, `podGC`, template arguments; matches no repo version | Retire candidate (`bed-split-template` is the current one) |

**For the wave-2 bpview wiring:**
- `somatic-rna-fastq-v2`, `somatic-fusion-fastq-gpu-v2`, `somatic-rna-bam-v2` and `extract-reads-for-fusion-inspector-template` are live-ahead with small, consistent changes.
- Committing their live specs first lets the wiring PR start from what runs.
- `somatic-fusion-fastq-v2` needs a decision before it is wired.

## Live-only templates imported here (28)
`templates/` holds them verbatim: the live `spec`, with `status`, `uid`, `resourceVersion`, `managedFields`, `creationTimestamp` and the last-applied annotation removed. **None ran in the last two months.**

| Template | Created | Called by | Calls |
|---|---|---|---|
| `arriba-fusion-detection` | 2026-02-10 | `fusion-detection-ffpe` | — |
| `cnv-calling` | 2026-02-10 | — | `annotsv-template`, `cnvkit-analysis`, `cnvkit-annotate`, `cnvkit-export`, `cnvkit-filter`, `collect-somatic-data-templates`, `duckdb-template`, `filter-annotsv-template`, `oncokb-templates`, `scarhrd` |
| `cnv-fusion-combined` | 2026-02-11 | — | — |
| `cnvkit-analysis` | 2026-02-10 | `cnv-calling` | — |
| `cnvkit-annotate` | 2026-02-10 | `cnv-calling` | — |
| `cnvkit-export` | 2026-02-10 | `cnv-calling` | — |
| `cnvkit-filter` | 2026-02-10 | `cnv-calling` | — |
| `collect-qc-data-templates` | 2023-11-23 | `germline-rna-gpu-workflow`, `msi-tmb-coverage-workflow-template`, `somatic-testing-rna-gpu-workflow`, `somatic-testing-workflow`, `tps-workflow`, `wes-workflow`, `wgs-workflow` | — |
| `filter-and-annot-workflow` | 2024-07-16 | — | `bcftools-template`, `cravat-templates` |
| `fusion-detection-ffpe` | 2026-02-10 | — | `arriba-fusion-detection`, `fusion-report-generation`, `oncokb-fusion-annotation`, `starfusion-detection` |
| `fusion-report-generation` | 2026-02-10 | `fusion-detection-ffpe` | — |
| `get-node-rel-dirs` | 2023-06-16 | — | — |
| `get-rel-dirs` | 2023-06-19 | — | — |
| `gridss-template` | 2025-03-18 | — | — |
| `manta-calls` | 2024-01-28 | — | `annotsv-template`, `bcftools-template`, `duckdb-template`, `elprep-template`, `filter-annotsv-template`, `manta-template`, `oncokb-templates`, `samtools-template` |
| `no-filter-mutect2-calls` | 2024-01-27 | — | `annovar-template`, `annovar-to-maf-template`, `bcftools-template`, `cancervar-template`, `duckdb-template`, `mutect2-template`, `oncokb-templates`, `vlod-template` |
| `oncokb-fusion-annotation` | 2026-02-10 | `fusion-detection-ffpe` | — |
| `parabricks-bam2fastq-template` | 2024-07-20 | — | — |
| `parabricks-fastq2bam-template` | 2024-07-19 | — | — |
| `parabricks-template` | 2024-07-14 | — | — |
| `pgs-workflow` | 2024-12-15 | — | `check-file-path`, `copy-inputs-to-scratch`, `copy-temp-to-output`, `dir-clean-up`, `gennext-prs-template`, `get-abs-pgs-ref-dirs`, `update-genotype-snps-template` |
| `recompute-cds-coverage-workflow-template` | 2026-02-12 | — | `bed-validator-template`, `collect-qc-data-v2-templates`, `mosdepth-template`, `multiqc-template` |
| `scarhrd` | 2026-02-11 | `cnv-calling` | — |
| `starfusion-detection` | 2026-02-10 | `fusion-detection-ffpe` | — |
| `structural-variant-analysis` | 2025-03-19 | — | — |
| `survivor` | 2025-03-18 | — | — |
| `tiddit` | 2025-03-18 | — | — |
| `update-run-status-template` | 2023-10-28 | — | — |

Notes:
- **`collect-qc-data-templates`** is called by seven DAGs that are on `main` (`wes-workflow`, `wgs-workflow`, `tps-workflow`, …), so `main` was not self-contained without it.
- **`fusion-detection-ffpe` family:** `arriba-fusion-detection`, `starfusion-detection`, `fusion-report-generation` and `oncokb-fusion-annotation` form one pipeline. `arriba-fusion-detection` still disables Arriba's `duplicates` filter (#40, item 5).
- **`cnv-calling` / `cnvkit-*` / `scarhrd`** form another pipeline. **`cnv-fusion-combined`** submits `cnv-calling` and `fusion-detection-ffpe` as child workflows.
- **Four no longer pass `argo lint` against today's live templates**, so they cannot run as they are:
  - `pgs-workflow` calls `get-abs-pgs-ref-dirs`, which exists neither on the cluster nor in the repo.
  - `manta-calls` doesn't pass `sample-id`, which `manta-template` now requires.
  - `no-filter-mutect2-calls` calls `vlod-template.vlod-merge-with-vcf`, which no longer exists.
  - `filter-and-annot-workflow` doesn't pass `regions-file` to `bcftools-template`.

  The other 24 pass.
- **Node names:** `cnv-calling`, `cnv-fusion-combined` and `fusion-detection-ffpe` hard-code a node name as a parameter default, as `cnv-workflow-template` and `msi-tmb-coverage-workflow-template` on `main` already do. They are kept verbatim here. Making them parameters without a default is one of the proposed fixes.

Proposal for each group: keep the two pipelines, `collect-qc-data-templates` and `recompute-cds-coverage-workflow-template` if their owners still need them. Retire the rest, starting with the four that no longer lint, after an owner check.

## Live-only templates not imported (23)
- **Tests and examples (13), retire candidates:** `complete-workflow-test`, `example-template-3ffd`, `example-template-3fv`, `fan-out-in-params-workflow`, `fastqc-v2-template` (used only by `test-workflow-v2`), `hostpath-example`, `multiqc-template-test`, `sprig-os-base-test`, `subtract-test`, `test-workflow`, `test-workflow-v2`, `trim-galore-test`, `verifybamid-template-test`. Their specs stay on record in `snapshots/live-2026-11-13/`.
- **Another project (9 templates):** out of scope, not imported.
- **`caller-workflow`, a retire candidate:** not imported. It is unused, and its parameter defaults carry test sample data.

## Repo-only templates (15)
These are on `main` but not on the cluster, and no live template calls them: `arriba-to-oncokb-template`, `batch-orchestrator-prototype`, `calc-tmb-template`, `caller-dv-workflow`, `caller-hc-workflow`, `filter-manta-template`, `freebayes-template`, `multicaller-template`, `parabricks-germline-pipeline-template`, `parallel-caller-workflow`, `per-sample-somatic-fusion`, `regenerate-qc-workflow`, `sample-worker-prototype`, `somatic-testing-dna-workflow`, `wrapper-batch-somatic-fusion`. Proposal: keep the prototypes that are in progress; move the rest to a retired folder or delete them, after an owner check.

## Images without a pinned tag (live)
The two named in #33 are `trinityctat/starfusion:latest` (`fusion-inspector-template`, 147 runs) and `qldrsc/duckdb:latest` (`duckdb-template`, 487 runs). Others in use:

| Image | Templates in use (runs) |
|---|---|
| `eeacms/rsync:latest` | `copy-temp-to-output` (533) |
| `griffan/verifybamid2:latest` | `verifybamid-template` (533) |
| `blairy/requests` (no tag) | `change-sample-state-template` (533), `post-process-cnv` (9), `post-process-variant-counter` (6) |
| `alpine:latest` / `alpine` | `dir-clean-up` (542), `copy-inputs-to-scratch` (62), `remove-temp-fqs-template` (62) |
| `ghcr.io/seqeralabs/rustqc:latest` | `rustqc-template` (147) |
| `alpine/kubectl` (no tag) | `get-node-config` (9), `get-base-node-config` (4) |

The Parabricks templates take their image from `{{workflow.parameters.parabricksImage}}`, so the pin belongs to the caller.

## Proposed follow-up PRs (not in this one)
1. **Commit live → repo** for the live-ahead templates in use: `copy-temp-to-output`, `parabricks-bqsr-template`, `extract-reads-for-fusion-inspector-template`, `somatic-rna-fastq-v2`, `somatic-rna-bam-v2`, `somatic-fusion-fastq-gpu-v2`, `copy-inputs-to-scratch`, `exomedepth-template`, `post-process-cnv`. Nothing is applied; afterwards these match live exactly. This unblocks the wave-2 bpview wiring.
2. **Apply repo → live** (gate:live-template): `filter-msisensor-pro-template`, `trim-galore-template`, `annotsv-to-oncokb-template`, and `fastqc-template` with its workflow-level fields dropped.
3. **Pin images by digest**, starting with `fusion-inspector` and `duckdb`. `duckdb-template` waits for the OncoKB track's change to that file.
4. **Owner decisions:** `somatic-fusion-fastq-v2`, `msi-tmb-coverage-workflow-template`, `wrapper-workflow-prs`, the imported live-only templates, the retire candidates and the repo-only templates.
5. **Token in pod logs:** the live `update-run-status-template` and `collect-qc-data-templates` print the cluster token in their pod logs, and `collect-qc-data-templates` posts to the gennext-dev backend. Stop printing the token, and point it at the right backend or retire it.
6. **Node-name defaults:** make them required parameters in `cnv-workflow-template`, `msi-tmb-coverage-workflow-template` and the imported templates that have them.
