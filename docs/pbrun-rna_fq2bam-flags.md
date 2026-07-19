# pbrun rna_fq2bam — Complete Flag List (Parabricks 4.3.1-1)

## Input/Output File Options

| Flag | Default | Description |
|------|---------|-------------|
| `--ref REF` | — | Path to the reference file. |
| `--in-fq [IN_FQ ...]` | — | Paired-end FASTQ files, optionally followed by a quoted read group string (e.g. `@RG\tID:foo\tLB:lib1\tPL:bar\tSM:sample\tPU:foo`). Repeatable. |
| `--in-se-fq [IN_SE_FQ ...]` | — | Single-ended FASTQ files, optionally followed by a quoted read group string. Repeatable. |
| `--genome-lib-dir GENOME_LIB_DIR` | — | Path to a STAR genome resource library directory (must be pre-indexed by user). |
| `--output-dir OUTPUT_DIR` | — | Directory that will hold all generated files. |
| `--out-bam OUT_BAM` | — | Path of the output BAM file. |

## Tool Options

### Output / Read Groups

| Flag | Default | Description |
|------|---------|-------------|
| `--out-prefix OUT_PREFIX` | — | Prefix filename for output data. |
| `--read-files-command READ_FILES_COMMAND` | — | Command to decompress input (e.g. `zcat`, `bzcat`). |
| `--read-group-sm READ_GROUP_SM` | — | SM (sample) tag for read groups. |
| `--read-group-lb READ_GROUP_LB` | — | LB (library) tag for read groups. |
| `--read-group-pl READ_GROUP_PL` | — | PL (platform) tag for read groups. |
| `--read-group-id-prefix READ_GROUP_ID_PREFIX` | — | Prefix used for ID and PU tags; an identifier is appended per FASTQ pair. |

### STAR Alignment

| Flag | Default | Description |
|------|---------|-------------|
| `--num-sa-bases NUM_SA_BASES` | 14 | Length of SA pre-indexing string (10–15 recommended). For small genomes scale down to `min(14, log2(GenomeLength)/2 - 1)`. |
| `--max-intron-size MAX_INTRON_SIZE` | 0 | Max intron size; 0 → `(2^winBinNbits)*winAnchorDistNbins`. |
| `--min-intron-size MIN_INTRON_SIZE` | 21 | Min intron size; smaller gaps treated as deletions. |
| `--min-match-filter MIN_MATCH_FILTER` | 0 | Min matched bases required for alignment output. |
| `--min-match-filter-normalized MIN_MATCH_FILTER_NORMALIZED` | 0.66 | Same as above, normalized to read length. |
| `--out-filter-intron-motifs OUT_FILTER_INTRON_MOTIFS` | None | `None`, `RemoveNoncanonical`, or `RemoveNoncanonicalUnannotated`. |
| `--max-out-filter-mismatch MAX_OUT_FILTER_MISMATCH` | 10 | Max mismatches allowed for an alignment to be output. |
| `--max-out-filter-mismatch-ratio MAX_OUT_FILTER_MISMATCH_RATIO` | 0.3 | Max ratio of mismatches to mapped length. |
| `--max-out-filter-multimap MAX_OUT_FILTER_MULTIMAP` | 10 | Max loci a read may map to before all alignments are dropped. |
| `--out-reads-unmapped OUT_READS_UNMAPPED` | None | `None` or `Fastx` (writes `Unmapped.out.mate1/2`). |
| `--out-sam-unmapped OUT_SAM_UNMAPPED` | None | `None`, `Within`, `KeepPairs`, or `Within_KeepPairs`. |
| `--out-sam-attributes OUT_SAM_ATTRIBUTES [...]` | Standard | Any combination of `{NH, HI, AS, nM, NM, MD, jM, jI, XS, MC, ch}`, or `None` / `Standard` / `All`. |
| `--out-sam-strand-field OUT_SAM_STRAND_FIELD` | None | `None` or `intronMotif` (Cufflinks-like strand flag). |
| `--out-sam-mode OUT_SAM_MODE` | Full | `None`, `Full`, or `NoQS` (no quality scores). |
| `--out-sam-mapq-unique OUT_SAM_MAPQ_UNIQUE` | 255 | MAPQ value for unique mappers (0–255). |
| `--min-score-filter MIN_SCORE_FILTER` | 0.66 | Min alignment score, normalized to read length. |
| `--min-spliced-mate-length MIN_SPLICED_MATE_LENGTH` | 0.66 | Min mapped length for a spliced mate, normalized to mate length. |
| `--max-junction-mismatches M M M M` | [0, -1, 0, 0] | Max mismatches for splice-junction stitching: (1) non-canonical, (2) GT/AG & CT/AC, (3) GC/AG & CT/GC, (4) AT/AC & GT/AT. -1 = no limit. |
| `--max-out-read-size MAX_OUT_READ_SIZE` | 100000 | Max SAM record size per read (bytes). |
| `--max-alignments-per-read MAX_ALIGNMENTS_PER_READ` | 10000 | Max alignments per read considered. |
| `--score-gap SCORE_GAP` | 0 | Splice junction penalty (motif-independent). |
| `--seed-search-start SEED_SEARCH_START` | 50 | Search start point through the read; split pieces will not exceed this length. |
| `--max-bam-sort-memory MAX_BAM_SORT_MEMORY` | 0 | Max RAM (bytes) for BAM sorting; 0 → use genome index size. |
| `--align-ends-type ALIGN_ENDS_TYPE` | Local | `Local` (soft-clip allowed) or `EndToEnd` (no soft-clipping). |
| `--align-insertion-flush ALIGN_INSERTION_FLUSH` | None | `None` or `Right` (flush ambiguous insertions right). |
| `--max-align-mates-gap MAX_ALIGN_MATES_GAP` | 0 | Max gap between mates; 0 → derived from `(2^winBinNbits)*winAnchorDistNbins`. |
| `--min-align-spliced-mate-map MIN_ALIGN_SPLICED_MATE_MAP` | 0 | Min mapped length for a spliced mate. |
| `--max-collapsed-junctions MAX_COLLAPSED_JUNCTIONS` | 1000000 | Max number of collapsed junctions. |
| `--min-align-sj-overhang MIN_ALIGN_SJ_OVERHANG` | 5 | Min overhang for spliced alignments. |
| `--min-align-sjdb-overhang MIN_ALIGN_SJDB_OVERHANG` | 3 | Min overhang for annotated (sjdb) spliced alignments. |
| `--sjdb-overhang SJDB_OVERHANG` | 100 | Donor/acceptor sequence length on each side of junctions (ideally `mate_length − 1`). |

### Chimeric Detection

| Flag | Default | Description |
|------|---------|-------------|
| `--min-chim-overhang MIN_CHIM_OVERHANG` | 20 | Min overhang for `Chimeric.out.junction`. |
| `--min-chim-segment MIN_CHIM_SEGMENT` | 0 | Min chimeric segment length; 0 disables chimeric output. |
| `--max-chim-multimap MAX_CHIM_MULTIMAP` | 0 | Max chimeric multi-alignments; 0 → legacy unique-only scheme. |
| `--chim-multimap-score-range CHIM_MULTIMAP_SCORE_RANGE` | 1 | Score range for multi-mapping chimeras (only with `--max-chim-multimap > 1`). |
| `--chim-score-non-gtag CHIM_SCORE_NON_GTAG` | -1 | Penalty for non-GT/AG chimeric junction. |
| `--min-non-chim-score-drop MIN_NON_CHIM_SCORE_DROP` | 20 | Threshold drop in best non-chimeric score (relative to read length) to trigger chimeric detection. |
| `--out-chim-format OUT_CHIM_FORMAT` | 0 | 0 = no headers; 1 = trailing comment lines (command line, Nreads). |
| `--out-chim-type OUT_CHIM_TYPE` | None | `Junctions`, `WithinBAM`, `WithinBAM_HardClip`, or `WithinBAM_SoftClip`. |
| `--two-pass-mode TWO_PASS_MODE` | None | `None` or `Basic` (insert first-pass junctions on the fly). |

### Duplicate Marking

| Flag | Default | Description |
|------|---------|-------------|
| `--no-markdups` | off | Skip Mark Duplicates step; return BAM after sorting. |
| `--read-name-separator SEP [SEP ...]` | `/` | Character(s) separating the trimmed portion of read names. |

## Performance Options

| Flag | Default | Description |
|------|---------|-------------|
| `--num-threads NUM_THREADS` | 4 | Worker threads per GPU. |

## Run Options

| Flag | Default | Description |
|------|---------|-------------|
| `-h, --help` | — | Show help message and exit. |
| `--verbose` | off | Enable verbose output. |
| `--x3` | off | Show full command-line arguments. |
| `--logfile LOGFILE` | stderr only | Path to log file. |
| `--tmp-dir TMP_DIR` | `.` | Directory for temporary files. |
| `--with-petagene-dir WITH_PETAGENE_DIR` | — | PetaGene install path; requires PetaLink preloaded via `LD_PRELOAD`. Honors env vars `PETASUITE_REFPATH`, `PGCLOUD_CREDPATH`, `PetaLinkMode` (set to `+write` to emit compressed BAM/FASTQ). |
| `--keep-tmp` | off | Retain temporary files after completion. |
| `--no-seccomp-override` | off | Do not override Docker seccomp options. |
| `--version` | — | Print compatible software versions. |
| `--num-gpus NUM_GPUS` | 1 | Number of GPUs to use. |
