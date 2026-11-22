# Fusion testing ground (prototype)

**Prototype, not part of any production pipeline.** Nothing calls these templates and they carry no
orchestrator label. They are applied and submitted by hand, with the maintainer's go, on an idle
node, for test runs in the DEV database only. Part of #29 / #30.

## `fusions-proto-pg`

Builds the new fusion event data (docs/fusions/event-contract.md, contract v1) for one past RNA
sample, from the outputs that sample already has, without touching them:

1. `check-inputs`: refuses non-UUID run or sample IDs, unusual sample names, and a dev run folder that
   overlaps the source sample's folder.
2. `copy-legacy`: copies the source's `results/fusions/{arriba,star-fusion,fusion-report,
   fusion-inspector,oncokb}`, `results/db/{s}_fusions.duckdb` and `bam/{s}_RNA.dedup.bam(.bai)` into each
   dev sample's folder, renaming the sample prefix, so the legacy fusion view works for the dev
   sample. It runs once per dev sample (a marker file), so a retry never overwrites a loaded DuckDB.
   Arm b's BAM is a hard link to arm a's copy.
3. Per Arriba arm, one dev sample each:
   - `a`: `max_sensitivity` (today's settings);
   - `b`: `filters_restored`, today's thresholds with Arriba's `blacklist`, `low_entropy` and
     `multimappers` filters on (Q5). Pass `.` for dev sample b to run arm a only.

   Arriba (own step: the live template can't express arm b) → `fusion-report` and
   `filter-artifact-fusions` (the live templates) → `fusion-events rescue` → `bpview` →
   `fusion-events build` → `load-events`. STAR-Fusion and OncoKB results are reused from the source:
   STAR's chimeric junctions are no longer kept, and OncoKB is annotated per gene pair.
4. `load-events` writes `fusion_events` and `fusion_events_meta` into the dev sample's copied fusions
   DuckDB in one transaction (meta last). Nothing else in that file changes.
5. Sample state (`running`, `completed`, `error`) and, at the end, run state go to the **dev**
   backend, `gennext-backend.gennext-dev`. The address is fixed in the template; there is no
   parameter for it.

Outputs per dev sample (under its output folder):

| Path | What |
|---|---|
| `results/fusions/bpview/{s}_bpview.*` | contigs, reads, exons, events.json |
| `results/fusions/fusion-events/{s}_fusion_events.tsv`, `_meta.tsv` | the event table as loaded |
| `results/fusions/proto/` | the arm's Arriba, fusion-report and filter outputs, the rescued Arriba TSV, step timings |
| `results/db/{s}_fusions.duckdb` | the copy, plus `fusion_events` and `fusion_events_meta` |

Images are pinned by digest. `fusion-report-template` and `filter-artifact-fusions-template` are
used as they are applied in the cluster.

## Running it

Gated: applying the template, creating the dev run (`POST /api/v0/runs`, never `start`) after a dev
database backup, and submitting the workflow each need the maintainer's go. The submit parameters
are `node`, `source-sample-dir`, `source-sample-id`, `dev-user-id`, `dev-run-id`,
`dev-sample-id-a`/`-uuid-a` and optionally `dev-sample-id-b`/`-uuid-b`. Sample IDs stay out of this
repository.

## Removing a test

Delete the dev run in the dev backend, then the dev run's output folder on the node, and the
workflow. The source sample is only read.
