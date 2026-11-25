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
| Legacy (unchanged) | `results/fusions/fusion-inspector/{s}_finspector.*`, `results/fusions/fusion-report/{s}_fusions.json` (or unprefixed `fusions.json`), table `advanced_combined_fusions` | today's pipeline |

`fusion_events_meta` is a two-column table (`key VARCHAR, value VARCHAR`):

| key | value |
|---|---|
| `contract_version` | `1` |
| `builder_version` | fusion-events version, e.g. `0.4.0` |
| `prototype` | `true` while the rules are the prototype's |
| `arriba_settings` | `max_sensitivity` or `filters_restored` (Q5); see §5 |
| `pon` | `none`, or the PoN file name and size |
| `oncokb_exact_pair` | the rule that decided `oncokb_pair_specific` (it states the rule; no data file is involved) |
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
| `legacy` | not `events`, and any of: the table `advanced_combined_fusions`; `results/fusions/fusion-inspector/{s}_finspector.fa.fai` (only that path, so a stray file elsewhere can't change the mode); the fusion-report JSON that today's table route (`GET /fusions`) reads, `results/fusions/fusion-report/{s}_fusions.json` or the unprefixed `results/fusions/fusion-report/fusions.json`, resolved the same way | today's table and FusionInspector view |
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
| `event_rank` | INTEGER | Default order, 1 = first: tier; then, from fusion-events 0.6.0, events that are a curated database pair (`curated_known`: Mitelman, COSMIC) or a pair OncoKB answers for specifically (`oncokb_pair_specific`), at a `canonical_junction`, with at least 1 unique molecule; then supported before `low_support`, unique molecules (descending), caller priority. |
| `sample_id` | VARCHAR | Sample name |
| `fusion` | VARCHAR | `GENE5--GENE3`, the join key to the gene-pair data |
| `gene5`, `gene3` | VARCHAR | Display genes |
| `gene5_aliases`, `gene3_aliases` | VARCHAR[] | Every name the callers gave each side |
| `breakpoint5`, `breakpoint3` | VARCHAR | `chr:pos:strand` (1-based). When `orientation_known` is false, which side is 5′ and the strands are Arriba's guess from its read directions. |
| `orientation_known` | BOOLEAN | false: an Arriba call without strand information |

### Structure (new in v1)
| Field | Type | Meaning |
|---|---|---|
| `transcript5`, `transcript3` | VARCHAR | The MANE Select transcript used for exon numbers |
| `exon5` | INTEGER | Last exon of the 5′ partner kept in the fusion (number in `transcript5`) |
| `exon3` | INTEGER | First exon of the 3′ partner kept in the fusion |
| `site5`, `site3` | VARCHAR | Where the breakpoint lies in that transcript: `exon_boundary`, `exon`, `intron` or `outside` |
| `reading_frame` | VARCHAR | Arriba's: `in-frame`, `out-of-frame`, `stop-codon`. When no caller gives one (rescued Arriba calls, STAR-Fusion-only events), derived from the MANE Select CDS phase (fusion-events ≥ 0.6.0), only when both breakpoints are at exon boundaries inside the CDS; never `stop-codon` (stop codons at the junction aren't detected). Null when not computable |
| `reading_frame_source` | VARCHAR | `arriba` or `mane` (derived; the UI says so); null when `reading_frame` is null. Tiers use the caller's frame only |
| `canonical_junction` | BOOLEAN | Both breakpoints at exon boundaries of the MANE Select transcripts (`site5` = `site3` = `exon_boundary`) |
| `adjacent_readthrough` | BOOLEAN | Same chromosome and strand, partners under 50 kb apart (`READTHROUGH_MAX_DISTANCE`), and no other gene (with a MANE Select transcript) overlapping the interval between the breakpoints |
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

**Rule: every support count shown matches the reads view.** Every count the view, the detail panel
or the report shows as support (`unique_molecules`, `molecules_junction`, `molecules_spanning`,
`pcr_copies`, `start_sites`) equals what can be counted in the event's contig (`bpview_contig`) in
the bpview BAM that the read view and desktop IGV load. The unit is the fragment: one QNAME.
- **Fragment class:** the BAM tag `XG`, the same on both mates. `molecules_junction` = non-duplicate
  fragments with `XG = 1junction`; `molecules_spanning` = non-duplicate fragments with `XG = 2spanning`;
  `unique_molecules` = their sum.
- **PCR copy:** a fragment where any record has the duplicate flag (0x400). `pcr_copies` = duplicate
  `1junction` + `2spanning` fragments. Duplicates of other classes don't count.
- **`start_sites`:** distinct `reference_start` values of the non-duplicate `1junction` fragments,
  taking per fragment the leftmost record whose own class (`XF`) is `junction`.
- **Shown, not counted:**
  - The evidence classes `junction_weak` (`XG = 3junction_weak`) and `chimeric_other` (`6chimeric_other`)
    are drawn in full, copies included, but are not support. Their numbers are in events.json
    (`fragments_nondup`, `fragments_dup`).
  - The context classes `wt5`, `wt3` and `other` are drawn up to `max_context` (40) non-duplicate
    fragments each, without copies. events.json counts all of them, before the cap. So a context count
    is not a drawn count, and the view should show context reads without a number or label the count
    "counted, up to 40 shown".
- **Events bpview didn't draw** (`bpview_contig` null; only past bpview's event cap or for a contig
  outside the genome FASTA) keep the callers' counts (`support_source` `arriba` or `starfusion`). They
  are labelled as caller counts and have no reads view, so the rule doesn't apply to them.
- **Caller numbers** (`arriba_*`, `starfusion_*`) stay caller columns, labelled as such.
- **Per-caller counts in the BAM** (bpview ≥ 0.4.0, events.json `bpview_version`): `XC:Z` on both mates
  of every fragment a caller cites (Arriba's `read_identifiers`, STAR-Fusion's `JunctionReads` and
  `SpanningFrags`) is `arriba`, `starfusion` or `arriba,starfusion`, always in that order; a fragment no
  caller cites has no `XC`. Cited fragments are always drawn, copies included and past `max_context`.
  - `{caller}_bam_molecules`: fragments with `XC` containing the caller, `XG` ∈ {`1junction`,
    `2spanning`}, and no record with 0x400.
  - `{caller}_bam_pcr_copies`: the same with at least one record flagged 0x400.
  - **Junction reads (bpview ≥ 0.4.1)** are either a chimeric split read (two records, primary +
    supplementary) or one spliced record whose blocks lie on both partners, joined by an N/D gap from one
    breakpoint to the other: STAR's alignment of a read-through between neighbouring genes. Both are drawn
    crossing the junction and counted the same way. Earlier bpview drew the spliced kind as wild-type and
    didn't count it.
  - `{caller}_bam_split_reads_5p` / `_3p`: the non-duplicate `1junction` fragments among those, by the
    side of the split read. `XJ:Z` (`5p` or `3p`) on a junction-class record (`XF = junction`) is the
    partner of its primary, non-supplementary alignment: Arriba's "anchor". For a spliced read-through
    record, which has no supplementary part, it's the partner holding more of its aligned bases. A fragment's side is its
    first junction mate's (read 1's `XJ`, else read 2's). Spanning pairs have no side, so 5p + 3p can be
    below the molecules.
  - **Caveat:** a cited fragment that bpview classes differently (e.g. `3junction_weak`, or a context
    class) carries `XC` but isn't counted. The view shows it as "cited, not counted". Both callers' name
    lists include duplicates; the molecule fields count non-duplicates only, and the duplicate-inclusive
    numbers are in events.json.
  - events.json, per event: `callers_bam` = `{arriba|starfusion: {fragments_nondup: {class: n},
    fragments_dup: {class: n}, split_5p, split_3p, cited, cited_missing}}`. `cited` is the number of
    distinct names the caller lists; `cited_missing` counts those not among the event's fragments (outside
    its segments, or no mate projected). They are never in the BAM.
  - `arriba_split_reads_5p` / `_3p` are Arriba's **own** `split_reads1`/`2` mapped to 5′/3′: caller numbers,
    not countable in the BAM.
  - A caller that cites no read names has null `{caller}_bam_*` fields (unknown, not zero). This is the
    case for every **rescued** Arriba call: Arriba writes `read_identifiers` only for kept calls.
  - Arriba's own numbers and its `_bam_` numbers can differ for kept calls too. Arriba cites reads that
    bpview classes as wild-type, for example a soft-clipped read whose clipped part has no supplementary
    alignment, and it counts by its own rules. Show both, each with its source.
- **IGV display:** desktop IGV and igv.js hide duplicates by default and downsample deep regions.
  The view has to switch downsampling off, or raise its depth, for its counts to match. The BAM has
  no QC-fail, secondary or supplementary records. Some junction reads have MAPQ 0, which IGV shows faded.

| Field | Type | Meaning |
|---|---|---|
| `unique_molecules` | INTEGER | Fragments supporting the junction, counted as `support_source` says |
| `molecules_junction`, `molecules_spanning` | INTEGER | The two parts of `unique_molecules` (bpview only) |
| `start_sites` | INTEGER | Distinct start positions of the non-duplicate junction reads (bpview only); for evidence thresholds such as ≥3 start sites (Q13) |
| `pcr_copies` | INTEGER | Duplicate junction and spanning fragments (bpview only) |
| `support_source` | VARCHAR | `bpview`: bpview's non-duplicate fragments, the same count whichever caller found the event. Otherwise `arriba`: Arriba's split reads + discordant mates after its own dedup. Otherwise `starfusion`: STAR-Fusion's junction + spanning reads, which are **not necessarily deduplicated**; label them "reads", not "molecules". |
| `low_support` | BOOLEAN | `unique_molecules` below the minimum (default 2). Shown, never dropped; capped at tier 2 (§5) |

### Annotation
| Field | Type | Meaning |
|---|---|---|
| `known_db` | VARCHAR[] | Databases recording the pair (fusion-report's, plus Arriba's Mitelman tag) |
| `curated_known` | BOOLEAN | In Mitelman or COSMIC |
| `oncokb_oncogenic`, `oncokb_mutation_effect` | VARCHAR | OncoKB, for the gene pair |
| `oncokb_highest_level`, `oncokb_highest_dx_level` | VARCHAR | e.g. `LEVEL_1`, `LEVEL_Dx1` |
| `oncokb_pair_specific` | BOOLEAN | **OncoKB's answer for this specific gene pair**: OncoKB's annotation is about this pair, not only about one gene's fusions, and not a loss-of-function match. Decided by a rule over the OncoKB API response the annotator stored for this fusion at annotation time, never from a copy of OncoKB's data. The rule: if the response records the alterations OncoKB matched (`ONCOKB_MATCHED_ALTERATIONS`, stored by the annotator from gennext-workflows#75, e.g. `EML4-ALK Fusion;Fusions`) and it isn't empty, the pair is specific when one entry is exactly "GENE5-GENE3 Fusion" or "GENE3-GENE5 Fusion". That answer decides on its own, even for a loss-of-function gene, because it is then explicitly about this pair. Otherwise (outputs from before #75, or OncoKB named nothing): `VARIANT_IN_ONCOKB` = True, two different genes, `MUTATION_EFFECT` not Loss-of-function / Likely Loss-of-function, `ONCOGENIC` not Unknown. The effect test is needed because OncoKB also sets `VARIANT_IN_ONCOKB` for tumour-suppressor and intragenic fusions it matches as truncating. Null when the pair has no OncoKB annotation (e.g. a rescued pair fusion-report never saw) |
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
| `tier_flags` | VARCHAR[] | What held it back: `out_of_frame`, `weak_caller_evidence`, `low_unique_molecules`, `low_support`, `rescued`, `rescue_capped`, `no_molecules`, `readthrough_capped`, `orientation_unknown`, `pon_recurrent`, `pon_capped` |
| `tier_reason` | VARCHAR | The same in words, for the hover and the detail |

### Caller and annotator details (fusion-events 0.5.0, additive)
Raw values for the detail panel, appended after the columns above; `contract_version` stays 1. They are
**caller and annotator values**, labelled as such, never the event's read count (see "Support"). Null
when that caller didn't report the event; lists `[]` when empty.

**Arriba**, from the arm's own Arriba row (the kept or rescued call). Sides 1/2 are Arriba's
`gene1`/`gene2` as written, which are the reverse of the event's 5′/3′ when an unoriented call was
flipped, so show them with `arriba_gene1`/`arriba_gene2`.

| Field | Type | Arriba column |
|---|---|---|
| `arriba_gene1`, `arriba_gene2` | VARCHAR | `gene1`, `gene2` (between genes: `A(dist),B(dist)`) |
| `arriba_breakpoint1`, `arriba_breakpoint2` | VARCHAR | `breakpoint1`, `breakpoint2` |
| `arriba_strand1`, `arriba_strand2` | VARCHAR | `strand1(gene/fusion)`, `strand2(gene/fusion)` |
| `arriba_split_reads1`, `arriba_split_reads2` | INTEGER | `split_reads1`, `split_reads2` (their sum is `arriba_split_reads`) |
| `arriba_coverage1`, `arriba_coverage2` | INTEGER | `coverage1`, `coverage2`: Arriba's coverage, a caller number; label it as Arriba's |
| `arriba_site1`, `arriba_site2` | VARCHAR | `site1`, `site2` |
| `arriba_filters_raw` | VARCHAR | `filters` as written, e.g. `duplicates(448)` |
| `arriba_duplicates` | INTEGER | the n of `duplicates(n)`: **Arriba's own** duplicate count (0 when absent). It differs from `pcr_copies` (bpview's count from the BAM); label it "Arriba's duplicate count", never next to "PCR copies" without its source |
| `arriba_tags` | VARCHAR[] | `tags` |
| `arriba_retained_domains1`, `arriba_retained_domains2` | VARCHAR[] | `retained_protein_domains`, split at `\|` and `,` |
| `arriba_closest_genomic_breakpoint1`, `…2` | VARCHAR | `closest_genomic_breakpoint1`, `…2` |
| `arriba_gene_id1`, `arriba_gene_id2`, `arriba_transcript_id1`, `arriba_transcript_id2` | VARCHAR | Arriba's genes and chosen transcripts (may differ from `transcript5`/`transcript3`, MANE) |
| `arriba_fusion_transcript`, `arriba_peptide_sequence` | VARCHAR | `fusion_transcript`, `peptide_sequence` (`\|` marks the junction); capped at 2,000 characters, then `...[truncated, N characters]` |

**STAR-Fusion**, from the event's STAR-Fusion row (left = 5′, right = 3′; for a junction reported
twice, the row kept). `starfusion_junction_reads`, `starfusion_spanning_frags` and `starfusion_ffpm`
are above.

| Field | Type | STAR-Fusion column |
|---|---|---|
| `starfusion_left_gene`, `starfusion_right_gene` | VARCHAR | `LeftGene`, `RightGene` (e.g. `FGFR3^ENSG…`) |
| `starfusion_left_breakpoint`, `starfusion_right_breakpoint` | VARCHAR | `LeftBreakpoint`, `RightBreakpoint` |
| `starfusion_est_junction`, `starfusion_est_spanning` | DOUBLE | `est_J`, `est_S`: STAR-Fusion's multimapping-weighted estimates, caller numbers; label them as STAR-Fusion's |
| `starfusion_splice_type` | VARCHAR | `SpliceType` (raw; `splice_type` above is the event's) |
| `starfusion_large_anchor_support` | BOOLEAN | `LargeAnchorSupport`: `YES_LDAS` true, `NO_LDAS` false |
| `starfusion_left_break_dinuc`, `starfusion_right_break_dinuc` | VARCHAR | `LeftBreakDinuc`, `RightBreakDinuc` |
| `starfusion_left_break_entropy`, `starfusion_right_break_entropy` | DOUBLE | `LeftBreakEntropy`, `RightBreakEntropy` |
| `starfusion_annots` | VARCHAR[] | `annots`, parsed as api-edge parses them |

**OncoKB** and **fusion-report**, per gene pair from `fusions_with_oncokb.json`:

| Field | Type | Source |
|---|---|---|
| `oncokb_highest_sensitive_level`, `oncokb_highest_resistance_level`, `oncokb_highest_px_level` | VARCHAR | `HIGHEST_SENSITIVE_LEVEL`, `HIGHEST_RESISTANCE_LEVEL`, `HIGHEST_PX_LEVEL` |
| `oncokb_treatments` | VARCHAR[] | one entry per non-empty `LEVEL_*` (1–4, R1–R2, Dx1–3, Px1–3): `"LEVEL_3B: Erdafitinib"` |
| `oncokb_mutation_effect_citations`, `oncokb_treatment_citations`, `oncokb_dx_citations`, `oncokb_px_citations` | VARCHAR[] | `MUTATION_EFFECT_CITATIONS`, `TX_CITATIONS`, `DX_CITATIONS`, `PX_CITATIONS`, split at `;` (PMIDs, or abstract text with a URL) |
| `oncokb_gene_in_oncokb`, `oncokb_variant_in_oncokb` | BOOLEAN | `GENE_IN_ONCOKB`, `VARIANT_IN_ONCOKB` |
| `oncokb_variant_summary` | VARCHAR | `ONCOKB_VARIANT_SUMMARY`: OncoKB's variantSummary sentence (annotator ≥ gennext-workflows#75; null before) |
| `oncokb_matched_alterations` | VARCHAR[] | `ONCOKB_MATCHED_ALTERATIONS`, split at `;`: the alterations OncoKB matched, e.g. `EML4-ALK Fusion`, `Fusions` (empty before #75 or when OncoKB named none) |
| `pair_fii_explained` | VARCHAR | fusion-report's `Explained FII` |

**Both callers, side by side on the event's 5′/3′** (fusion-events 0.5.0 with bpview 0.4.0; INTEGER; null
when the caller didn't call the event, or for `*_bam_*`, when bpview didn't draw it). Definitions in
"Support" above.

| Field | Source |
|---|---|
| `arriba_split_reads_5p`, `arriba_split_reads_3p` | Arriba's own `split_reads1`/`2`, swapped when the call was flipped |
| `arriba_bam_split_reads_5p`, `arriba_bam_split_reads_3p`, `arriba_bam_molecules`, `arriba_bam_pcr_copies` | bpview BAM, fragments with `XC` ∋ `arriba` |
| `starfusion_bam_split_reads_5p`, `starfusion_bam_split_reads_3p`, `starfusion_bam_molecules`, `starfusion_bam_pcr_copies` | bpview BAM, fragments with `XC` ∋ `starfusion` (STAR-Fusion reports no per-side numbers itself) |
| `arriba_cited_missing`, `starfusion_cited_missing` | names the caller cites that bpview didn't find (events.json `callers_bam.*.cited_missing`) |

Not available from these inputs: OncoKB summary text (not in the annotator's output), and per-database
identifiers or counts (fusion-report gives database names only; they are `known_db`). The OncoKB JSON's
own `arriba` block comes from the run that made the JSON and isn't used.

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
  The database and OncoKB hits come from the sample's `fusions_with_oncokb.json`, which annotates only
  the pairs fusion-report saw (Arriba's kept calls and STAR-Fusion's). For a discarded-only pair that
  branch therefore rarely applies; Arriba's Mitelman tag (from its known-fusions file) and the
  STAR-Fusion branch do the work.
- **Rescue rule (Berk, 25 Nov; fusion-events 0.6.0).** A rescued call stays above tier 3 only when OncoKB
  gives its answer for that specific gene pair (`oncokb_pair_specific`, a rule over the OncoKB API answer; see §4), whatever its molecule count; partner-agnostic
  levels or a database tag alone don't lift it. Even then it never reaches tier 1. `tier_flags` gets
  `rescue_capped`, and `tier_reason` says which part applied.
- **Adjacent-gene read-through (Berk, 25 Nov; fusion-events 0.6.0).** An `adjacent_readthrough` event is
  tier 3 unless OncoKB's answer is for that specific pair (`oncokb_pair_specific`; `tier_flags` `readthrough_capped`; `tier_reason`
  gives the distance and the 50 kb limit). Intrachromosomal fusions further apart, or with a gene between
  them (e.g. GOPC–ROS1, ~240 kb), are untouched.
- **No counted molecule (Berk, 25 Nov; fusion-events 0.6.0).** An event with 0 `unique_molecules` is tier 3,
  whatever the path, a kept curated call included (`tier_flags` `no_molecules`).
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
- `contract_version` changes only for a breaking change. api-edge treats a version it doesn't know
  as if the sample had no event data (§2): `legacy` when the legacy files exist, otherwise `none`.
