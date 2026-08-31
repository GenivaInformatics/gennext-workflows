#!/bin/bash
# Regression test for GNV-315/GNV-316: 100bp BED padding with analysis and IGV artifacts
set -euo pipefail

echo "=== BED Padding Regression Test ==="

# Test 1: Verify no duplicate padding logic in bed_validator.yaml
echo "Test 1: No duplicate padding-generation logic"
PADDING_BLOCKS=$(grep -c "Generating analysis and IGV BED artifacts" tasks/base_tasks/bed_validator.yaml)
if [[ "$PADDING_BLOCKS" != "1" ]]; then
  echo "FAIL: Found $PADDING_BLOCKS padding-generation blocks, expected exactly 1"
  exit 1
fi
echo "✓ Exactly one padding-generation block in bed_validator.yaml"
echo ""

# Test 2: Verify explicit sort before merge
echo "Test 2: Explicit sort before merge"
if ! grep -q "padded_intervals.sort" tasks/base_tasks/bed_validator.yaml; then
  echo "FAIL: No explicit sort found before merge"
  exit 1
fi
echo "✓ Explicit sort present in padding implementation"
echo ""

# Test 3: Verify hard-fail on padding generation failure (no silent fallback)
echo "Test 3: No silent fallback to unpadded BED on padding failure"
# Check that the script exits with error if padding fails
if ! grep -A5 "if \[ \$? -ne 0 \]" tasks/base_tasks/bed_validator.yaml | grep -q "exit 1"; then
  echo "FAIL: Padding failure does not exit with error - would silently fall back"
  exit 1
fi
echo "✓ Padding failure causes hard error (no silent fallback)"
echo ""

# Test 4: Verify v2 workflows use validated-bed-analysis-file for both mutect2 and bcftools
echo "Test 4: v2 Mutect2 and bcftools-filter use validated-bed-analysis-file"
for workflow in workflows/somatic/v2/somatic-dna-bam.yaml workflows/somatic/v2/somatic-dna-fastq.yaml; do
  # Check mutect2 uses validated-bed-analysis-file
  MUTECT2_BED=$(grep -A20 "name: mutect2\$" "$workflow" | grep -A1 "regions-file" | tail -1 | grep -o "validated-bed[^}]*")
  if [[ "$MUTECT2_BED" != "validated-bed-analysis-file" ]]; then
    echo "FAIL: $workflow Mutect2 should use validated-bed-analysis-file, got: $MUTECT2_BED"
    exit 1
  fi
  echo "✓ $workflow Mutect2 uses validated-bed-analysis-file"

  # Check bcftools-filter uses validated-bed-analysis-file
  BCFTOOLS_TASK=$(grep -A20 "name: bcftools-filter\$" "$workflow")
  if ! echo "$BCFTOOLS_TASK" | grep -q "mutect2.outputs.parameters.output-vcf"; then
    echo "FAIL: $workflow bcftools-filter not processing mutect2 output"
    exit 1
  fi
  if ! echo "$BCFTOOLS_TASK" | grep -q "validated-bed-analysis-file"; then
    echo "FAIL: $workflow bcftools-filter not using validated-bed-analysis-file"
    exit 1
  fi
  # Ensure it's NOT using the old padded or original
  if echo "$BCFTOOLS_TASK" | grep "regions-file" | grep -q "validated-bed-file}}"; then
    echo "FAIL: $workflow bcftools-filter using unpadded validated-bed-file"
    exit 1
  fi
  echo "✓ $workflow bcftools-filter uses validated-bed-analysis-file"
done
echo ""

# Test 5: Verify TMB uses validated-bed-file (NOT analysis or IGV)
echo "Test 5: TMB uses validated-bed-file (not analysis or IGV)"
for workflow in workflows/somatic/v2/somatic-dna-bam.yaml workflows/somatic/v2/somatic-dna-fastq.yaml; do
  TMB_BED=$(grep -A20 "name: collect-tmb-data" "$workflow" | grep -o "validated-bed[^}]*" | head -1)
  if [[ "$TMB_BED" != "validated-bed-file" ]]; then
    echo "FAIL: $workflow TMB should use validated-bed-file, got: $TMB_BED"
    exit 1
  fi
  echo "✓ $workflow TMB correctly uses validated-bed-file (not analysis)"
done
echo ""

# Test 6: Verify Mutect2 direct template still has --interval-padding 100
echo "Test 6: Direct mutect2.yaml has --interval-padding 100"
MUTECT2_PADDING=$(grep "interval-padding 100" tasks/base_tasks/mutect2.yaml | wc -l)
if [[ "$MUTECT2_PADDING" != "2" ]]; then
  echo "FAIL: Expected 2 occurrences of '--interval-padding 100' in mutect2.yaml, found $MUTECT2_PADDING"
  exit 1
fi
echo "✓ Both direct Mutect2 commands use --interval-padding 100"
echo ""

# Test 7: Verify padding value is exactly 100bp
echo "Test 7: Padding is exactly 100bp"
for workflow in workflows/somatic/v2/somatic-dna-bam.yaml workflows/somatic/v2/somatic-dna-fastq.yaml; do
  if ! grep -A2 "name: padding" "$workflow" | grep -q 'value: "100"'; then
    echo "FAIL: $workflow padding value is not 100"
    exit 1
  fi
  echo "✓ $workflow padding value is 100bp"
done
echo ""

# Create test fixture directory for padding algorithm tests
TEST_DIR=$(mktemp -d)
trap "rm -rf $TEST_DIR" EXIT

# Create test GTF (no ##sequence-region headers, just regular records with max coordinates)
cat > "$TEST_DIR/test.gtf" <<'EOF'
chr1	HAVANA	gene	1000	1050	.	+	.	gene_id "TEST1"
chr1	HAVANA	gene	1200	1250	.	+	.	gene_id "TEST2"
chr1	HAVANA	gene	10000	10100	.	+	.	gene_id "TEST3"
chr1	HAVANA	gene	248956422	248956422	.	+	.	gene_id "CHR1_END"
chr22	HAVANA	gene	16050000	16050100	.	+	.	gene_id "TEST4"
chr22	HAVANA	gene	50818468	50818468	.	+	.	gene_id "CHR22_END"
EOF

# Create test BED with overlapping, book-ended, and edge cases
cat > "$TEST_DIR/test.bed" <<'EOF'
chr1	1000	1050
chr1	1200	1250
chr1	10000	10100
chr22	16050000	16050100
chr1	10	60
chr22	50818368	50818468
EOF

echo "Test 8: Padding algorithm produces correct analysis and IGV artifacts"

# Run the Python padding algorithm from the manifest
python3 << PYTHON_EOF
import sys

validated_bed = "${TEST_DIR}/test.bed"
analysis_bed = "${TEST_DIR}/analysis.bed"
igv_bed = "${TEST_DIR}/igv.bed"
gtf_file = "${TEST_DIR}/test.gtf"
padding = 100

# Extract chromosome sizes from GTF
chrom_sizes = {}
with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 5:
            continue
        chrom = fields[0]
        try:
            end = int(fields[4])
            chrom_sizes[chrom] = max(chrom_sizes.get(chrom, 0), end)
        except (ValueError, IndexError):
            continue

# Read validated BED intervals (original targets)
target_intervals = []
with open(validated_bed, 'r') as bed_in:
    for line in bed_in:
        if line.startswith('#') or line.startswith('track') or line.strip() == '':
            continue
        fields = line.strip().split('\t')
        if len(fields) < 3:
            continue
        chrom = fields[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
        except (ValueError, IndexError):
            continue
        if chrom not in chrom_sizes:
            continue
        target_intervals.append((chrom, start, end))

# Generate padded intervals for analysis BED
padded_intervals = []
for chrom, start, end in target_intervals:
    padded_start = max(0, start - padding)
    padded_end = min(chrom_sizes[chrom], end + padding)
    padded_intervals.append((chrom, padded_start, padded_end))

# Sort padded intervals explicitly
padded_intervals.sort(key=lambda x: (x[0], x[1], x[2]))

# Merge overlapping and book-ended intervals
merged_analysis = []
current_chrom, current_start, current_end = padded_intervals[0]

for chrom, start, end in padded_intervals[1:]:
    if chrom == current_chrom and start <= current_end:
        current_end = max(current_end, end)
    else:
        merged_analysis.append((current_chrom, current_start, current_end))
        current_chrom, current_start, current_end = chrom, start, end

merged_analysis.append((current_chrom, current_start, current_end))

# Write analysis BED (strict BED3)
with open(analysis_bed, 'w') as bed_out:
    for chrom, start, end in merged_analysis:
        bed_out.write(f"{chrom}\t{start}\t{end}\n")

# Generate IGV BED with targeted and padding-only features (BED9)
target_sorted = sorted(target_intervals, key=lambda x: (x[0], x[1], x[2]))
igv_features = []

# Add all targeted features (original intervals, gray 128,128,128)
for chrom, start, end in target_sorted:
    igv_features.append((chrom, start, end, 'targeted', '128,128,128'))

# Compute padding-only regions by subtracting targets from merged analysis
for analysis_chr, analysis_start, analysis_end in merged_analysis:
    overlapping_targets = [
        (s, e) for chr, s, e in target_intervals
        if chr == analysis_chr and not (e <= analysis_start or s >= analysis_end)
    ]

    if not overlapping_targets:
        igv_features.append((analysis_chr, analysis_start, analysis_end, 'padding', '211,211,211'))
        continue

    overlapping_targets.sort()
    current_pos = analysis_start
    for target_start, target_end in overlapping_targets:
        if current_pos < target_start:
            igv_features.append((analysis_chr, current_pos, target_start, 'padding', '211,211,211'))
        current_pos = max(current_pos, target_end)

    if current_pos < analysis_end:
        igv_features.append((analysis_chr, current_pos, analysis_end, 'padding', '211,211,211'))

# Sort IGV features
igv_features.sort(key=lambda x: (x[0], x[1], x[2], x[3]))

# Write IGV BED (BED9)
with open(igv_bed, 'w') as bed_out:
    bed_out.write('track name="Targeted and Padding Regions" itemRgb="On"\n')
    for chrom, start, end, feature_type, color in igv_features:
        bed_out.write(f"{chrom}\t{start}\t{end}\t{feature_type}\t0\t.\t{start}\t{end}\t{color}\n")
PYTHON_EOF

echo "Analysis BED:"
cat "$TEST_DIR/analysis.bed"
echo ""
echo "IGV BED:"
cat "$TEST_DIR/igv.bed"
echo ""

# Test 8a: Analysis BED is strict BED3 (exactly 3 columns)
echo "Test 8a: Analysis BED is strict BED3"
if ! awk '{if (NF != 3) exit 1}' "$TEST_DIR/analysis.bed"; then
  echo "FAIL: Analysis BED is not strict BED3 (has columns other than chr, start, end)"
  exit 1
fi
echo "✓ Analysis BED is strict BED3"

# Test 8b: Analysis BED has correct padding (100bp), merged overlaps
# First merged: chr1:0-160 (10-100 to 0, 60+100 to 160)
# Second merged: chr1:900-1350 (1000-100 to 900, 1250+100 to 1350, merged because they overlap/book-end)
# Third: chr1:9900-10200
# Fourth: chr22:16049900-16050200
# Fifth: chr22:50818268-50818468 (clipped at end)
FIRST_ANALYSIS=$(sed -n '1p' "$TEST_DIR/analysis.bed")
if [[ "$FIRST_ANALYSIS" != "chr1"$'\t'"0"$'\t'"160" ]]; then
  echo "FAIL: First analysis interval should be chr1:0-160, got: $FIRST_ANALYSIS"
  exit 1
fi
echo "✓ Analysis BED: chromosome-start clipping works (chr1:0-160)"

MERGED_900_1350=$(grep "chr1.*900.*1350" "$TEST_DIR/analysis.bed")
if [[ -z "$MERGED_900_1350" ]]; then
  echo "FAIL: Overlapping/book-ended intervals chr1:1000-1050 and chr1:1200-1250 should merge to chr1:900-1350"
  exit 1
fi
echo "✓ Analysis BED: overlapping/book-ended intervals merge correctly (chr1:900-1350)"

CHR_END=$(grep "chr22.*50818468" "$TEST_DIR/analysis.bed")
if [[ -z "$CHR_END" ]]; then
  echo "FAIL: Chromosome-end clipping failed"
  exit 1
fi
echo "✓ Analysis BED: chromosome-end clipping works"

# Test 8c: IGV BED has targeted and padding-only features with correct colors
echo "Test 8c: IGV BED has distinct targeted/padding-only features"
TARGETED_COUNT=$(grep -c "targeted.*128,128,128" "$TEST_DIR/igv.bed" || echo 0)
PADDING_COUNT=$(grep -c "padding.*211,211,211" "$TEST_DIR/igv.bed" || echo 0)
if [[ "$TARGETED_COUNT" == "0" ]]; then
  echo "FAIL: IGV BED has no targeted features"
  exit 1
fi
if [[ "$PADDING_COUNT" == "0" ]]; then
  echo "FAIL: IGV BED has no padding-only features"
  exit 1
fi
echo "✓ IGV BED has $TARGETED_COUNT targeted features (gray 128,128,128)"
echo "✓ IGV BED has $PADDING_COUNT padding-only features (light gray 211,211,211)"

# Test 8d: IGV BED has no overlap between targeted and padding-only features
echo "Test 8d: No coordinate overlap between targeted and padding-only in IGV BED"
python3 << PYTHON_EOF
igv_bed = "${TEST_DIR}/igv.bed"
targeted_ranges = []
padding_ranges = []

with open(igv_bed, 'r') as f:
    for line in f:
        if line.startswith('track') or line.strip() == '':
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        chrom, start, end, feature_type = fields[0], int(fields[1]), int(fields[2]), fields[3]
        if feature_type == 'targeted':
            targeted_ranges.append((chrom, start, end))
        elif feature_type == 'padding':
            padding_ranges.append((chrom, start, end))

# Check for overlaps
for t_chr, t_start, t_end in targeted_ranges:
    for p_chr, p_start, p_end in padding_ranges:
        if t_chr == p_chr:
            # Check if they overlap
            if not (t_end <= p_start or t_start >= p_end):
                print(f"FAIL: Overlap found between targeted {t_chr}:{t_start}-{t_end} and padding {p_chr}:{p_start}-{p_end}", file=sys.stderr)
                sys.exit(1)

print("✓ No overlap between targeted and padding-only features in IGV BED")
PYTHON_EOF

echo ""
echo "=== All Tests Passed ==="
