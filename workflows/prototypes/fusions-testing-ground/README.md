# Fusion testing ground (prototype)

**Prototype, not part of any production pipeline.** Nothing calls these templates and they carry no
orchestrator label. They are applied and submitted by hand, with the maintainer's go, on an idle
node, for test runs in the DEV database only. Part of #29 / #30.

## `fusions-proto-pg`

Builds the new fusion event data (docs/fusions/event-contract.md, contract v1) for one past RNA
sample, from the outputs that sample already has, without touching them:

1. `check-inputs`: refuses non-UUID run or sample IDs, unusual sample names, and a dev run folder that
   overlaps the source sample's folder.
2. `copy-legacy` (only after the dev backend accepted both samples as part of the dev run): copies the source's `results/fusions/{arriba,star-fusion,fusion-report,
   fusion-inspector,oncokb}`, `results/db/{s}_fusions.duckdb` and `bam/{s}_RNA.dedup.bam(.bai)` into each
   dev sample's folder, renaming the sample prefix, so the legacy fusion view works for the dev
   sample. It runs once per dev sample: a marker file records the source, so a retry never
   overwrites a loaded DuckDB, and a resubmit with another source is refused. It refuses a target
   folder that already has content without that marker (a real sample's).
   Arm b's BAM is a hard link to arm a's copy.
3. Per Arriba arm, one dev sample each:
   - `a`: `max_sensitivity` (today's settings);
   - `b`: `filters_restored`, today's thresholds with Arriba's `blacklist`, `low_entropy` and
     `multimappers` filters on (Q5). Pass `.` for dev sample b to run arm a only.

   Arriba (own step: the live template can't express arm b) → `fusion-report` and
   `filter-artifact-fusions` (the live templates) → `fusion-events rescue` → `bpview` →
   `fusion-events build` → `load-events`. bpview draws every merged event (no candidate list), so
   every count the event table reports can be counted in the BAM the read view and IGV load. STAR-Fusion and OncoKB results are reused from the source:
   STAR's chimeric junctions are no longer kept, and OncoKB is annotated per gene pair.
4. `load-events` writes `fusion_events` and `fusion_events_meta` into the dev sample's copied fusions
   DuckDB in one transaction (meta last). Nothing else in that file changes.
5. Sample state (`running`; `completed` or `error` per arm) and, after a successful run, run state go to the **dev**
   backend, `gennext-backend.gennext-dev`. The address is fixed in the template; there is no
   parameter for it.

Outputs per dev sample (under its output folder):

| Path | What |
|---|---|
| `results/fusions/bpview/{s}_bpview.*` | contigs, reads, exons, events.json |
| `results/fusions/fusion-events/{s}_fusion_events.tsv`, `_meta.tsv` | the event table as loaded |
| `results/fusions/proto/` | the arm's Arriba, fusion-report and filter outputs, the rescued Arriba TSV, step timings |
| `results/db/{s}_fusions.duckdb` | the copy, plus `fusion_events` and `fusion_events_meta` |

This template's own steps pin their images by digest and set memory requests and limits. The two
reused live templates are the exception, as they are applied in the cluster today:
`fusion-report-template` (`quay.io/biocontainers/fusion-report:4.0.1--py313hdfd78af_0`, tag only, no
memory limit) and `filter-artifact-fusions-template` (`python:3.12-slim`, tag only, no memory limit).
Pinning them is #33's live → repo reconcile, not this prototype's.

## Running it

Gated: applying the template, creating the dev run (`POST /api/v0/runs`, never `start`) after a dev
database backup, and submitting the workflow each need the maintainer's go. The submit parameters
are `node`, `source-sample-dir`, `source-sample-id`, `dev-user-id`, `dev-run-id`,
`dev-sample-id-a`/`-uuid-a` and optionally `dev-sample-id-b`/`-uuid-b`. Sample IDs stay out of this
repository.

## Removing a test

Delete the dev run in the dev backend, then the dev run's output folder on the node, and the
workflow. The source sample is only read.
