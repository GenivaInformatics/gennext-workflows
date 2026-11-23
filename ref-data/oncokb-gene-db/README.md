# oncokbGene: ANNOVAR on OncoKB's isoforms

An ANNOVAR gene database with OncoKB's own transcript for every curated gene, so that
the protein change we compute is numbered the way OncoKB numbers it.

- **Curated genes:** OncoKB's `grch38Isoform`, with its version from the Ensembl
  release OncoKB pins.
- **Every other gene:** Ensembl MANE Select. So the database is a full gene set, like
  refGene and ensGene.
- **Gene names:** OncoKB's `hugoSymbol`. Spaces and parentheses become `_` for ANNOVAR,
  so `CDKN2A (p14)` is `CDKN2A_p14`.

On the 23 Nov 2026 comparison against OncoKB's own genomic answers (v7.7, Ensembl 111),
the protein change on these isoforms equals OncoKB's for 100% (narrow panel) and 99.84%
(wide panel) of coding variants in curated genes.

## Files

| file | |
|---|---|
| `build.sh` | builds the database and its tarball |
| `build_genepred.py` | Ensembl GTF → 16-column genePred |
| `build_isoforms.py` | the isoform table for annovar-to-maf, plus the translation check |
| `install.sh` | installs a tarball into one ANNOVAR database directory |
| `copy-to-node.sh` | runs `install.sh` on a node, through a short-lived pod |
| `make_fixture.py` | regenerates `tests/fixtures/oncokb-gene-db` from a build |

A build `oncokb-<data version>_ensembl-<release>` holds four files. Each is named by its
version, for example `hg38_oncokbGene.oncokb-v7.7_ensembl-111.txt`:

| file | |
|---|---|
| `hg38_oncokbGene.<v>.txt` | the genePred ANNOVAR reads |
| `hg38_oncokbGeneMrna.<v>.fa` | the transcript sequences (ANNOVAR's `retrieve_seq_from_fasta.pl`) |
| `hg38_oncokbGene.<v>.isoforms.tsv` | one row per OncoKB isoform: exons, Ensembl protein, 3′ UTR |
| `hg38_oncokbGene.<v>.version` | inputs and their sha256, the ANNOVAR image, counts, known exceptions |

The tarball also has `SHA256SUMS`.

## How the pipeline uses it

- **annovar-template:** with `oncokb-gene-db: auto` (the default), it adds `oncokbGene`
  as the last gene protocol, but only where `annovardb-dir` has a complete install:
  the four links all pointing at one version's files. Links that disagree (a
  half-finished switch) mean the protocol is skipped, with a message. `off` never adds
  it; any other value fails the step. This adds the `Func/Gene/GeneDetail/ExonicFunc/AAChange.oncokbGene`
  columns before `Otherinfo`. It also copies the isoform table and the version next to
  the multianno, as `<sample>.oncokbGene.isoforms.tsv` and `<sample>.oncokbGene.version`.
  A node without the database runs exactly as before.
- **annovar-to-maf:** maftools-rs 0.2.1 or later. When the isoform table is there, rows
  of curated genes show the protein change on OncoKB's isoform where it differs from
  MANE, and the MAF gets `*_mane` and `oncokb_*` columns (see maftools-rs's README).
  refhub's OncoKB fallback reads the `oncokb_isoform*` columns.

## Build

Needs docker and python3, and runs under the build caps:

```bash
systemd-run --user --scope -p MemoryMax=24G -p MemorySwapMax=0 -p CPUQuota=1600% \
  ./build.sh oncokb-v7.7-allCuratedGenes.json /path/to/hg38.analysisSet.fa /path/to/work v7.7 111
```

- **allCuratedGenes:** OncoKB's `/api/v1/utils/allCuratedGenes`, saved for the data
  version OncoKB answers with.
- **Reference:** the FASTA the pipeline aligns to.
- **Downloads:** the Ensembl GTF and proteins come from ftp.ensembl.org, and must match
  the sha256 pinned in `build.sh` for release 111 (pass `GTF_SHA256` / `PEP_SHA256` for
  another release).
- **Output:** `work/dist/oncokbGene-<v>.tar.gz`. The script prints its sha256.
- **Reproducible:** a rebuild from the same inputs is byte-identical, whatever the
  builder's umask (the tarball stores fixed modes, owners and times).

`build_isoforms.py` fails the build if an OncoKB isoform is missing from the release, or
if a coding isoform doesn't translate to its Ensembl protein. Differences at residue 1
(non-ATG starts) and incomplete CDS are allowed and listed as notes. For v7.7 /
Ensembl 111 those are MYC, WT1 and RYBP.

**Known exception:** RYBP's Ensembl 111 model has an incomplete ORF. ANNOVAR skips it, so
RYBP variants get no `oncokbGene` protein change.

**When to rebuild:**
- OncoKB changes a curated gene's isoform: compare `allCuratedGenes` between data
  releases.
- OncoKB moves its Ensembl pin.

Install the new version next to the old one. The links switch, and each sample records
the version it used (`oncokb_gene_db` in the MAF).

## Install on a node

The ANNOVAR database directory is the one annovar-template mounts: `ref-paths`
`$.other.annovar` under the node's reference root.

```bash
./copy-to-node.sh --context <kube-context> <node> <annovardb-host-path> oncokbGene-oncokb-v7.7_ensembl-111.tar.gz <sha256>
```

`--context` is required, so the cluster is explicit. The path must be absolute and end
in `/annovardb`.

**What it does:**
- runs `install.sh` in an alpine pod on that node (namespace `gennext-dev` by default,
  256Mi limit), with only that directory mounted. The pod has no service-account
  token, a read-only root filesystem, and all capabilities dropped but
  `DAC_OVERRIDE` (to write a directory it may not own);
- checks the tarball's sha256 and its `SHA256SUMS`;
- moves the four versioned files in, read-only. Files already there are never
  overwritten;
- then points `hg38_oncokbGeneMrna.fa`, `hg38_oncokbGene.isoforms.tsv`,
  `hg38_oncokbGene.version` and, last, `hg38_oncokbGene.txt` at them as relative
  symlinks, each replaced by a rename.

**A sample annotating during the install** sees either no database or a complete one.

**An interrupted install** (a killed pod, Ctrl-C) can leave the links pointing at two
versions. annovar-template then skips the protocol on that node, and running the same
install again completes the switch: `install.sh` first removes the staging directory
and temp links an earlier run left, and checks at the end that all four links name
the new version.

**A switch mid-run:** if a sample's ANNOVAR step started on one version and the links
switched before it finished, the step fails rather than mixing versions. annovar-template
has no retryStrategy: rerun it with `argo retry`, which uses the new version.

**Undo:**

```bash
./copy-to-node.sh --context <kube-context> --remove <node> <annovardb-host-path>
```

This removes the links, `hg38_oncokbGene.txt` first, and keeps the versioned files. From
then on, samples skip the protocol again.

Running `install.sh` directly on a node works the same way:

```bash
install.sh <tarball> <annovardb-dir> <sha256>
```

## Tests

- `tests/test_annovar_oncokb_gene.sh`: annovar-template and annovar-to-maf in their
  production images, on the fixture. Without the database (or with links naming two
  versions, or `off`), their output is byte-identical to the checked-in baseline in
  `tests/fixtures/oncokb-gene-db/expected`, which the templates wrote before
  oncokbGene existed (`REGEN_GOLDEN=1` rewrites it from that commit). It also checks
  an unknown parameter value and a switch mid-run, and, with `CANCERVARDB` set (the
  mise task fetches CancerVar's data at a pinned commit), CancerVar's output.
- `tests/test_oncokb_gene_db_install.sh`: `install.sh` under busybox, including a
  re-run after an interrupted switch.
