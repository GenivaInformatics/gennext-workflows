#!/usr/bin/env bash
# Regenerates the synthetic reference files under ref/ (bgzip + tabix, in the
# htslib image). The paths and names are the real ones under
# <ref-base-dir>/hg38/other; the contents are a handful of made-up records
# that match variants.json. Output is committed; rerun only to change it.
set -euo pipefail
cd "$(dirname "$0")"
rm -rf ref
mkdir -p ref/cosmic ref/common_variants/gnomad4 ref/common_variants/turkish_variome

contigs() {
  printf '##contig=<ID=%s>\n' "$@"
}

{
  echo '##fileformat=VCFv4.2'
  echo '##INFO=<ID=TARGETED_SCREEN_SAMPLE_COUNT,Number=1,Type=Integer,Description="Samples">'
  contigs 3 12 13 17 22
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  # COSMIC uses bare contig names. Every row of variants.json except
  # chr17:7676060 A>T is a COSMIC hit.
  printf '3\t179203600\tCOSV1\tA\tG\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=12\n'
  printf '3\t179203650\tCOSV2\tG\tC\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=3\n'
  printf '12\t25245300\tCOSV3\tC\tT\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=40\n'
  printf '12\t25245310\tCOSV4\tG\tA\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=8\n'
  printf '12\t25245320\tCOSV5\tT\tG\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=2\n'
  printf '13\t32336300\tCOSV6\tC\tT\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=5\n'
  printf '17\t7676000\tCOSV7\tC\tT\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=30\n'
  printf '17\t7676010\tCOSV8\tG\tA\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=9\n'
  printf '17\t7676020\tCOSV9\tA\tG\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=4\n'
  printf '17\t7676030\tCOSV10\tT\tC\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=4\n'
  printf '17\t7676040\tCOSV11\tC\tA\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=4\n'
  printf '17\t7676050\tCOSV12\tG\tT\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=4\n'
  printf '17\t7676099\tCOSV13\tGC\tG\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=6\n'
  printf '22\t28734500\tCOSV14\tC\tT\t.\t.\tTARGETED_SCREEN_SAMPLE_COUNT=7\n'
} > ref/cosmic/Cosmic_CompleteTargetedScreensMutant_v103_GRCh38.vcf

{
  echo '##fileformat=VCFv4.2'
  echo '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">'
  contigs chr3 chr12 chr13 chr17 chr22
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'chr12\t25245300\t.\tC\tT\t.\tPASS\tAF=0.05\n'
  printf 'chr12\t25245320\t.\tT\tG\t.\tPASS\tAF=0.001\n'
  printf 'chr17\t7676099\t.\tGC\tG\t.\tPASS\tAF=0.2\n'
} > ref/common_variants/gnomad4/gnomad4_1pct.vcf

{
  echo '##fileformat=VCFv4.2'
  echo '##INFO=<ID=TR_AF,Number=A,Type=Float,Description="Turkish Variome allele frequency">'
  contigs chr3 chr12 chr13 chr17 chr22
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'chr12\t25245310\t.\tG\tA\t.\tPASS\tTR_AF=0.02\n'
} > ref/common_variants/turkish_variome/turkish_variome_1pct.vcf

{
  echo '##fileformat=VCFv4.2'
  echo '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">'
  contigs chr3 chr12 chr13 chr17 chr22
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  printf 'chr12\t25245300\t.\tC\tT\t.\tPASS\tAF=0.05\n'
  printf 'chr12\t25245320\t.\tT\tG\t.\tPASS\tAF=0.001\n'
} > ref/common_variants/gnomad4/gnomad4.1_exomes_pass_ac10.vcf

docker run --rm --network none --user "$(id -u):$(id -g)" -v "$PWD/ref:/ref" \
  quay.io/biocontainers/htslib:1.17--h81da01d_2 \
  sh -c 'for f in $(find /ref -name "*.vcf"); do bgzip -f "$f" && tabix -p vcf "$f.gz"; done'
find ref -type f | sort
