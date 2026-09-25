# Somatic Annotation and Filter Fixes

Outstanding changes from the Illumina concordance review of six 2026 samples
(run `14bb91b8`, node `gennext-02`), 2026-08-31 → 2026-09-07.

Full analysis: https://claude.ai/code/artifact/b8557cb9-428e-40a6-85fc-81dfcfaedb9f

## Context

Thirteen reported discrepancies against Illumina resolved to three defects plus
three variants genuinely absent from the reads. One defect (a stale ANNOVAR
gene database on `gennext-02`) is already fixed. The two below are not.

---

## TODO 1 — `tasks/base_tasks/mutect2-sharded.yaml`

### 1a. Raise the clustered-events threshold

`FilterMutectCalls` is invoked with no `--max-events-in-region`, so GATK's
default of **2** applies. A confirmed TP53 `p.Glu224Ter` in M-1438-26 — 191/269
alt reads, 73 % VAF, TLOD 425, ROQ 93, strand-balanced — was rejected because
two neighbouring 2–3 % FFPE C>A artifacts pushed `ECNT` to 4.

**Change** (~line 141):

```diff
             gatk FilterMutectCalls \
               -R "$REF" \
               -V "$UNFILTERED" \
               --stats "$UNFILTERED.stats" \
               --orientation-bias-artifact-priors "$PRIORS" \
+              --max-events-in-region 4 \
               -O "$OUT_VCF"
```

**Why 4** — measured across all six samples on variants whose *only* filter is
`clustered_events` (n=2433), scoring high-confidence as AF ≥ 10 % and ROQ ≥ 90:

| `--max-events-in-region` | released | high-conf | low-conf | marginal precision |
|---|---|---|---|---|
| 2 (default) | 0 | 0 | 0 | — |
| 3 | 216 | 107 | 109 | 49.5 % |
| **4** | **941** | **394** | **547** | **39.6 %** |
| 5 | 1426 | 560 | 866 | 34.2 % |
| 6 | 1782 | 679 | 1103 | 33.4 % |
| 9 | 2433 | 838 | 1595 | 17.6 % |

The TP53 target sits at `ECNT=4`, so 4 is the smallest value that recovers it.
Marginal precision is still ~40 % at 4 and degrades from there; 4 → 5 admits 485
more calls to gain 166 high-confidence ones.

FFPE deamination products are **not** released by this change — they carry
`orientation` as a second filter from the read-orientation model
(`LearnReadOrientationModel` + `--orientation-bias-artifact-priors`), which is
already configured and operates independently of the event count.

### 1b. Reduce artifact generation (test alongside 1a, not instead of it)

`--pileup-detection-snp-alt-threshold 0.02` promotes 2 % events into calls,
which is what inflates `ECNT` in the first place. Two of the four events in the
M-1438-26 cluster exist only because of it; without them `ECNT` would be 2 and
the TP53 variant would pass at GATK's untouched default.

Consider raising the threshold, or disabling SNV pileup detection for FFPE
inputs. Needs a sensitivity check on low-VAF true positives before adopting.

---

## TODO 2 — `tasks/base_tasks/annovar-to-maf.yaml`

### 2a. Normalise ANNOVAR's `.` before `annovarToMaf`

maftools 2.18.0 derives `Hugo_Symbol` from `Gene.refGene` — which is correct —
then destroys it:

```r
Variant_Classification = ifelse(is.na(ExonicFunc.refGene),
                                annovar_values[Func.refGene],        # branch we want
                                annovar_values[ExonicFunc.refGene])  # branch taken
Hugo_Symbol = ifelse(Variant_Classification == "IGR", "Unknown", Hugo_Symbol)
Hugo_Symbol = ifelse(is.na(Hugo_Symbol), "Unknown", Hugo_Symbol)
```

ANNOVAR runs with `--nastring .`, so non-exonic rows carry a literal dot. R
reads `"."` as a string, `is.na()` is FALSE, the wrong branch is taken,
`annovar_values["."]` returns `NA`, and `NA == "IGR"` evaluates to `NA` — not
FALSE — so `ifelse` propagates `NA` and the final sweep writes `"Unknown"`.

`annovar_values` already contains `splicing = "Splice_Site"`,
`upstream = "5'Flank"`, `intronic = "Intron"`, `UTR5`, `UTR3`. They are simply
never reached.

**Change** — preprocess to `/tmp` (not the hostPath output dir), then point
`annovarToMaf` at the normalised copy:

```sh
IN=/mnt/input/{{=sprig.osBase(inputs.parameters['annovar-out-f'])}}
awk -F'\t' -v OFS='\t' 'NR==1{print;next}
                        {if ($9==".") $9="";
                         if ($6=="intron") $6="intronic";
                         print}' "$IN" > /tmp/normalised.txt

Rscript -e "library(maftools); \
maf = annovarToMaf(annovar='/tmp/normalised.txt', sep='\t', \
refBuild='{{inputs.parameters.ref-version}}', MAFobj=TRUE, \
basename='{{inputs.parameters.sample-id}}')"
```

Column 9 is `ExonicFunc.refGene`, column 6 is `Func.refGene`. Confirm `awk`
exists in `bergun/r-maftools:v0.0.1` before merging.

Prefer this over changing ANNOVAR's global `--nastring`, which would rewrite
every column in the multianno and has a far wider blast radius.

**Verified on M-1037-26's real multianno:**

| | before | after |
|---|---|---|
| TERT `chr5:1295135` | `Unknown` / `NA` | **`TERT` / `5'Flank`** |
| rows with `Hugo_Symbol = "Unknown"` | 1994 | 232 |
| rows with a `Variant_Classification` | 893 | 2800 |
| Splice_Site recovered | 0 | 17 |
| 5'Flank (promoter) recovered | 0 | 11 |

1762 of 2887 rows regain their gene symbol. The residual 232 are 145 genuinely
intergenic (`IGR` → `Unknown` is correct maftools behaviour) and 87 carrying
ANNOVAR's `intron` label, which the `$6` substitution above clears.

### 2b. Fall back to `GeneDetail.refGene` for `tx` / `txChange`

Non-coding rows put HGVS in `GeneDetail.refGene`, not `AAChange.refGene`, so
`tx`/`exon`/`txChange`/`aaChange` stay empty even after 2a. Post-MANE the
fallback is unambiguous — one transcript per row.

### 2c. Synthesise promoter cDNA notation

ANNOVAR emits only `dist=N` for `upstream` variants, never a `c.` string. TERT
`C250T` is `c.-146C>T` because it is 67 bp upstream of the `NM_198253.3`
transcript start across a 79 bp 5'UTR. Either special-case the two hotspots
(`C228T` chr5:1,295,113 and `C250T` chr5:1,295,135) or derive `c.-N` from the
MANE GTF as `dist + 5'UTR length`.

---

## TODO 3 — node database parity

`gennext-02` was updated 2026-08-31 (copied from `acurare-01`, verified
byte-identical). Remaining:

| node | refGene | ensGene | version file | state |
|---|---|---|---|---|
| acurare-01, gennext-03, gennext-04 | 19342 `ca55c44418` | 19342 | 19342 | OK |
| gennext-02 | 19342 `ca55c44418` | 19342 | 19342 | fixed 08-31 |
| **gennext-01** | 19254 `c34df558e9` | **111108** | **76502** | partial |
| **ibg-02** | — | — | — | unverified (kubelet 503) |

Copy from `acurare-01`, not from the local `.mane` files — those are the
*unversioned* build (`11c3184f2179`, emits `NM_000546` not `NM_000546.6`).
`tar` preserves source ownership, so re-own `root:root 0644` after transfer or
the ANNOVAR container cannot read them.

---

## TODO 4 — guardrail

The MANE database makes `annovarToMaf`'s take-the-first-transcript behaviour
correct *by accident*: one transcript per gene means there is nothing else to
pick. If a comprehensive database is ever restored on any node the discrepancies
return silently.

Either pin a MANE Select lookup inside the conversion step, or assert at
pipeline start that `hg38_refGene.txt` carries one transcript per gene.

Cheap cross-check: `ONCOKB_PROTEIN_CHANGE` already carries canonical HGVS
independently (OncoKB queries by genomic coordinate). Surfacing it beside
`aaChange` would have made this visible in January.

---

## TODO 5 — re-issue affected reports

The delivered `.maf`, `_oncokb.maf` and `results.duckdb` for the six samples
still carry the old nomenclature; the database swap alone does not rewrite them.
Re-running `annovar` → `annovar-to-maf` → `oncokb` → `duckdb` regenerates them
with no re-alignment or re-calling needed.

M-1438-26 is the exception and the priority: TP53 `p.Glu224Ter` is a
likely-oncogenic LOF call at 73 % VAF currently absent from the report, and only
TODO 1a recovers it.

Also worth auditing how many other cases ran on `gennext-02` between 2026-01-26
and 2026-08-31 — every one carries the same non-canonical numbering.

**Note a side effect before re-issuing:** dropping to one transcript per gene
reclassifies 33–115 variants per sample from exonic to non-exonic (coding only
on a minor isoform — mostly `EGFR-AS1`, `PAX8-AS1`, `SETBP1`, with `BRCA1` and
`EGFR` entries becoming intronic). Correct behaviour for MANE-based reporting,
but coding-variant lists shrink and anyone diffing an old report against a new
one will see variants disappear.

---

## Not actionable — closed

- **RUNX1 `p.Ala357ProfsTer19`** (M-1652-26) — not in the reads. All eight
  `NM_001754` coding exons at 114–1244x; largest indel signal at any position is
  3 reads (~0.3 %). No Manta SV. Contamination 0.19 %.
- **RBM10 `c.967+2`** (M-1652-26) — all 23 donor and 23 acceptor sites covered
  829–1735x, none above 1 % alt. No RBM10 transcript has an exon boundary at
  c.967; RUNX1's exon-2 donor does (`chr21:34,799,299`, 260/260 reference).
- **TP53 `c.920-162`** (M-1642-26) — **coordinate cannot exist.** TP53 intron 8
  is 92 bp (`chr17:7,673,609-7,673,700`), so valid intronic positions are
  `c.919+1..c.919+46` and `c.920-46..c.920-1`. The only intron-8 signal is
  `c.920-1G>T` at chr17:7,673,609 (5 % VAF, `FILTER=orientation`) — almost
  certainly the intended variant, and it needs orthogonal confirmation rather
  than a pipeline change.
