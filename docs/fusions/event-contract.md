# Fusion event data contract (v1, prototype)

Status: **draft for the prototype** (`fusions-proto-pg`, F5 = api-edge#12, F6). Part of #29 / #30.
Readers: gennext-api-edge (F5) and gennext-frontend (F6). Producer: `fusion-events` in
`images/gennext-bpview` (0.4.0) and the `fusion_events` load step.

Rule for every reader: **a missing file, table, column or value means "unknown", never false or
zero.** Older samples have none of this. The prototype marks its own rules (`prototype_rules`), and
they may change after the scorecard.

## 1. Files and tables per sample

Paths are relative to the sample's output directory
(`{output-base}/{user}/{run}/{sample_uuid}/`, as api-edge resolves today). `{s}` is the sample name.

| What | Path | Written by |
|---|---|---|
| Event table | table `fusion_events` in `results/db/{s}_fusions.duckdb` | fusion_events load |
| Event metadata | table `fusion_events_meta` in the same DuckDB | fusion_events load |
| Event TSV (source of the table) | `results/fusions/fusion-events/{s}_fusion_events.tsv` | fusion-events build |
| bpview contigs | `results/fusions/bpview/{s}_bpview.fa`, `.fa.fai` | bpview |
| bpview reads | `results/fusions/bpview/{s}_bpview.bam`, `.bam.bai` | bpview |
| bpview exons | `results/fusions/bpview/{s}_bpview.bed.gz`, `.gtf.gz` (+ `.tbi`) | bpview |
| bpview events | `results/fusions/bpview/{s}_bpview.events.json` | bpview |
| Legacy (unchanged) | `results/fusions/fusion-inspector/{s}_finspector.*`, table `advanced_combined_fusions` | today's pipeline |

`fusion_events_meta` is a two-column table (`key VARCHAR, value VARCHAR`):

| key | value |
|---|---|
| `contract_version` | `1` |
| `builder_version` | fusion-events version, e.g. `0.4.0` |
| `prototype` | `true` while the rules are the prototype's |
| `arriba_settings` | `max_sensitivity` or `filters_restored` (Q5); see §5 |
| `pon` | `none`, or the PoN file name and size |
| `events` | number of rows in `fusion_events` |

The load step writes both tables in one transaction, `fusion_events_meta` last, so a sample never
has the meta table without a complete event table. The prototype loads them with its own step
(`fusions-proto-pg`); the F4 template `duckdb-fusion-events` loads an older column set and no
meta table, and is not used for contract v1.

## 2. Mode flag

Per sample, from what is on disk. api-edge computes it; nothing stores it.

| mode | When | Frontend |
|---|---|---|
| `events` | **both** tables exist: `fusion_events`, and `fusion_events_meta` with a `contract_version` api-edge knows | new event view |
| `legacy` | not `events`, and `advanced_combined_fusions` exists or `results/fusions/fusion-inspector/{s}_finspector.fa.fai` exists (only that path, so a stray file elsewhere can't change the mode) | today's table and FusionInspector view |
| `none` | neither | "no fusion data" state |

`bpview` (boolean) says separately whether the bpview files exist: all four of `{s}_bpview.fa`,
`.fa.fai`, `.bam` and `.bam.bai`. An `events` sample without
bpview files still gets the event view, but no reads view for any event.

## 3. api-edge routes (F5, additions only)

All under `/api/v1/variants/{run_id}/somatic`, with the same auth and the same `sampleName` /
`sampleUuid` query parameters as today's fusion routes. Today's `/fusions*` routes don't change.

| Route | Returns |
|---|---|
| `GET /fusions/mode` | `{ "mode": "events"\|"legacy"\|"none", "bpview": bool, "contractVersion": "1"\|null, "builderVersion": str\|null, "prototype": bool }` |
| `GET /fusions/events` | `{ "mode", "bpview", "contractVersion", "builderVersion", "prototype", "arribaSettings", "events": [event…] }`. A sample that is not in `events` mode returns 200 with `events: []` and its mode, never an error. |
| `GET /fusions/bpview/events` | the sample's `{s}_bpview.events.json`, as written |
| `GET /fusions/bpview/contigs` | `{ "contigs": [{ "name", "length", "eventId" }] }` from the `.fa.fai`, joined to the events |
| `GET /fusions/bpview/contig/fasta?contig=` | one contig's FASTA |
| `GET /fusions/bpview/contig/bam?contig=` | one contig's reads as a BAM (igv.js), as the FusionInspector route does |
| `GET /fusions/bpview/contig/bam/index?contig=` | its BAI |
| `GET /fusions/bpview/contig/bed?contig=` | exon BED records for the contig (tabix) |
| `GET /fusions/bpview/contig/gtf?contig=` | exon GTF records for the contig (tabix) |

An unknown contig returns 400 and a missing bpview file 404, as the FusionInspector routes do.
Contig names contain characters such as `-` and `:` and may contain others, so the caller URL-encodes
`contig` and api-edge decodes it (today's FusionInspector IGV path doesn't encode; the new one must).

Casing is deliberate: the response wrapper's own keys are camelCase like every api-edge response
(`contractVersion`, `arribaSettings`), while the event objects keep the DuckDB column names
(snake_case), so one field has one name in the table, the API and the filter body.
`events` are sorted in the default order (§4, `event_rank`).

## 4. Event fields

One object per breakpoint event. **Keys are the DuckDB column names, unchanged** (snake_case), so a
field has one name in the table, the API and the saved-filter body. List fields are JSON arrays
(`VARCHAR[]` in DuckDB), `[]` when empty, never null. Any other `null` = unknown.

Rows are written 5′ → 3′. One junction is one event: A–B and B–A never appear as two rows for
the same junction.

### Identity and order
| Field | Type | Meaning |
|---|---|---|
| `event_id` | VARCHAR | Stable key for selection, views and the report: the junction, `chr5:pos5:strand5>chr3:pos3:strand3` (strand `+`/`-`, or `.` on both sides when `orientation_known` is false). The same junction in a re-run gets the same ID. It is unique **within a sample** only (the same junction in two samples has the same ID), so selection and the report key on (`sample_uuid`, `event_id`). In the rare case of two events at one junction, the second gets `#2`. |
| `event_rank` | INTEGER | Default order, 1 = first: tier, then supported before `low_support`, then unique molecules (descending), then caller priority. |
| `sample_id` | VARCHAR | Sample name |
| `fusion` | VARCHAR | `GENE5--GENE3`, the join key to the gene-pair data |
| `gene5`, `gene3` | VARCHAR | Display genes |
| `gene5_aliases`, `gene3_aliases` | VARCHAR[] | Every name the callers gave each side |
| `breakpoint5`, `breakpoint3` | VARCHAR | `chr:pos:strand` (1-based) |
| `orientation_known` | BOOLEAN | false: an Arriba call without strand information |

### Structure (new in v1)
| Field | Type | Meaning |
|---|---|---|
| `transcript5`, `transcript3` | VARCHAR | The MANE Select transcript used for exon numbers |
| `exon5` | INTEGER | Last exon of the 5′ partner kept in the fusion (number in `transcript5`) |
| `exon3` | INTEGER | First exon of the 3′ partner kept in the fusion |
| `site5`, `site3` | VARCHAR | Where the breakpoint lies in that transcript: `exon_boundary`, `exon`, `intron` or `outside` |
| `reading_frame` | VARCHAR | Arriba's: `in-frame`, `out-of-frame`, `stop-codon`; null without Arriba |
| `splice_type` | VARCHAR | STAR-Fusion's `ONLY_REF_SPLICE` / `INCL_NON_REF_SPLICE`, or derived from Arriba's sites |
| `kinase_gene` | VARCHAR | The partner with a protein kinase domain (Pfam PF00069); the 3′ partner when both have one; null when neither has one |
| `kinase_retained` | BOOLEAN | At least 90% of the kinase partner's kinase domain lies on the kept side of its breakpoint (`kinase_domain_fraction`). Null when `kinase_gene` is null. |
| `kinase_domain_fraction` | DOUBLE | Fraction (0–1) of that kinase domain that is kept |

`exon5`/`exon3`/`site*` are null when the gene has no MANE Select transcript or the kept side runs
antisense to it. Kinase domains come from Arriba's bundled protein-domain annotation (Pfam),
the same file Arriba uses for `retained_protein_domains`.

### Callers
| Field | Type | Meaning |
|---|---|---|
| `callers` | VARCHAR[] | `arriba`, `starfusion` |
| `n_callers` | INTEGER | |
| `caller_arriba`, `caller_starfusion` | BOOLEAN | |
| `arriba_status` | VARCHAR | `kept` (Arriba reported it), `rescued` (brought back from Arriba's discarded list, §5) or null (no Arriba call) |
| `arriba_filters` | VARCHAR[] | For rescued calls: the Arriba filters that had removed it |
| `arriba_confidence` | VARCHAR | `high`, `medium`, `low` |
| `arriba_type` | VARCHAR | Arriba's event type |
| `arriba_split_reads`, `arriba_discordant_mates` | INTEGER | Arriba's counts (after its dedup) |
| `starfusion_junction_reads`, `starfusion_spanning_frags` | INTEGER | STAR-Fusion's counts |
| `starfusion_ffpm` | DOUBLE | |
| `breakpoint_offset` | INTEGER | With both callers: the larger distance (bp) between their breakpoints on either side; 0 = they agree exactly. Null with one caller. |

### Support
| Field | Type | Meaning |
|---|---|---|
| `unique_molecules` | INTEGER | Fragments supporting the junction, counted as `support_source` says |
| `molecules_junction`, `molecules_spanning` | INTEGER | The two parts of `unique_molecules` (bpview only) |
| `start_sites` | INTEGER | Distinct start positions of the non-duplicate junction reads (bpview only); for evidence thresholds such as ≥3 start sites (Q13) |
| `pcr_copies` | INTEGER | Duplicate fragments at the junction (bpview only) |
| `support_source` | VARCHAR | `bpview`: bpview's non-duplicate fragments, the same count whichever caller found the event. Otherwise `arriba`: Arriba's split reads + discordant mates after its own dedup. Otherwise `starfusion`: STAR-Fusion's junction + spanning reads, which are **not necessarily deduplicated**; label them "reads", not "molecules". |
| `low_support` | BOOLEAN | `unique_molecules` below the minimum (default 2). Shown, never dropped; capped at tier 2 (§5) |

### Annotation
| Field | Type | Meaning |
|---|---|---|
| `known_db` | VARCHAR[] | Databases recording the pair (fusion-report's, plus Arriba's Mitelman tag) |
| `curated_known` | BOOLEAN | In Mitelman or COSMIC |
| `oncokb_oncogenic`, `oncokb_mutation_effect` | VARCHAR | OncoKB, for the gene pair |
| `oncokb_highest_level`, `oncokb_highest_dx_level` | VARCHAR | e.g. `LEVEL_1`, `LEVEL_Dx1` |
| `pair_clinical_tier`, `pair_fii` | VARCHAR, DOUBLE | Today's gene-pair tier and FII, for comparison |

### Recurrence
| Field | Type | Meaning |
|---|---|---|
| `pon_gene_pair_count` | INTEGER | Other PoN samples with the same gene pair |
| `pon_junction_count` | INTEGER | Other PoN samples with the same junction (±10 bp) |
| `pon_size` | INTEGER | Other samples in the PoN |
| `pon_recurrent` | BOOLEAN | |

All null when the sample was built without a PoN (the prototype's default).

### Tier
| Field | Type | Meaning |
|---|---|---|
| `tier` | INTEGER | 1–3, a default prioritisation, never a filter. Every event is in the table. |
| `tier_basis` | VARCHAR | `oncokb_actionable`, `oncokb_oncogenic`, `curated_known`, `none` |
| `tier_flags` | VARCHAR[] | What held it back: `out_of_frame`, `weak_caller_evidence`, `low_unique_molecules`, `low_support`, `rescued`, `orientation_unknown`, `pon_recurrent`, `pon_capped` |
| `tier_reason` | VARCHAR | The same in words, for the hover and the detail |

### Links and provenance
| Field | Type | Meaning |
|---|---|---|
| `bpview_contig` | VARCHAR | The event's contig in the bpview files; null when bpview didn't draw it (then no reads view) |
| `builder_version` | VARCHAR | |
| `prototype_rules` | VARCHAR[] | Prototype rules that applied to this event: `rescue`, `low_support` |

The prototype's "selection", RNA QC and saved filters are not event fields; they stay with the
backend and the filter body.

## 5. Prototype rules (Q1, Q4, Q5)

- **Rescue (Q1).** A call from Arriba's `fusions.discarded.tsv` enters the table with
  `arriba_status = rescued` when all three hold: at least 1 unique read; every filter that removed it
  is `in_vitro`, `relative_support` or `min_support` (`duplicates(n)` is ignored); and the pair is known
  (Arriba's Mitelman tag, or a database or OncoKB hit for the pair) or STAR-Fusion calls the same
  junction within ±10 bp. Novel Arriba-only single-read calls stay out. bpview draws rescued events too.
- **Low support (Q1).** `low_support` when `unique_molecules` < 2. Capped at tier 2, ordered after
  supported events of the same tier. The view shows them by default; "Hide low support" is off by default.
- **Default view (Q4 A).** Event samples have no built-in caller filters. Sort by `event_rank`; tiers 1
  and 2 open, tier 3 collapsed.
- **Arriba arms (Q5).** The prototype runs Arriba twice on the same BAM and writes each arm to its own
  dev sample: `max_sensitivity` (today) and `filters_restored` (today's settings with Arriba's
  `blacklist`, `low_entropy` and `multimappers` filters back on). `fusion_events_meta.arriba_settings`
  says which.
- **Merge.** Arriba and STAR-Fusion calls are joined to the *nearest* call within ±10 bp, and STAR-Fusion
  rows at the same breakpoint are one event.

## 6. Compatibility

- Additions only: no existing DuckDB object, file or route changes. The legacy view keeps working
  for every sample, including the prototype's dev samples (their legacy files are copies).
- A newer contract adds fields; it never renames or retypes one. A reader ignores unknown fields.
- `contract_version` changes only for a breaking change, and api-edge then reports `legacy` for a
  version it doesn't know.
