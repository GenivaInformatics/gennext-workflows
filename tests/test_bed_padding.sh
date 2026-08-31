#!/bin/bash
# Regression test for GNV-315/GNV-316: 100bp BED padding
set -euo pipefail

# Use standard /bin/grep to avoid ugrep incompatibilities
#alias /bin/grep=/bin//bin/grep

echo "=== BED Padding Regression Test ==="

# Create test fixture directory
TEST_DIR=$(mktemp -d)
trap "rm -rf $TEST_DIR" EXIT

# Create a tiny test BED with overlapping/adjacent scenarios when padded
cat > "$TEST_DIR/test.bed" <<'EOF'
chr1	1000	1050
chr1	1200	1250
chr1	10000	10100
chr22	16050000	16050100
EOF

# Create a minimal .fai for chr1 and chr22
cat > "$TEST_DIR/test.fai" <<'EOF'
chr1	248956422	52	80	81
chr22	50818468	5136834921	80	81
EOF

echo "Input BED:"
cat "$TEST_DIR/test.bed"
echo ""

# Test 1: Verify bedtools slop + merge produces correct padding
echo "Test 1: 100bp padding with merge"
bedtools slop -i "$TEST_DIR/test.bed" -g "$TEST_DIR/test.fai" -b 100 | \
  bedtools merge -i - > "$TEST_DIR/padded.bed"

echo "Padded BED:"
cat "$TEST_DIR/padded.bed"
echo ""

# Verify first two intervals merged (1000-100=900, 1250+100=1350)
# Original: chr1:1000-1050 and chr1:1200-1250
# Padded: chr1:900-1150 and chr1:1100-1350 -> should merge to chr1:900-1350
FIRST_LINE=$(head -1 "$TEST_DIR/padded.bed")
if [[ "$FIRST_LINE" != "chr1"$'\t'"900"$'\t'"1350" ]]; then
  echo "FAIL: First two intervals should merge to chr1 900-1350"
  echo "Got: $FIRST_LINE"
  exit 1
fi

# Verify third interval padded correctly (10000-100=9900, 10100+100=10200)
SECOND_LINE=$(sed -n '2p' "$TEST_DIR/padded.bed")
if [[ "$SECOND_LINE" != "chr1"$'\t'"9900"$'\t'"10200" ]]; then
  echo "FAIL: Third interval should be chr1 9900-10200"
  echo "Got: $SECOND_LINE"
  exit 1
fi

# Verify fourth interval padded correctly and didn't merge
THIRD_LINE=$(sed -n '3p' "$TEST_DIR/padded.bed")
if [[ "$THIRD_LINE" != "chr22"$'\t'"16049900"$'\t'"16050200" ]]; then
  echo "FAIL: Fourth interval should be chr22 16049900-16050200"
  echo "Got: $THIRD_LINE"
  exit 1
fi

# Verify exactly 3 output intervals (2 merged, 1 separate, 1 separate)
LINE_COUNT=$(wc -l < "$TEST_DIR/padded.bed")
if [[ "$LINE_COUNT" != "3" ]]; then
  echo "FAIL: Expected 3 output intervals, got $LINE_COUNT"
  exit 1
fi

echo "✓ Padding produces exactly 100bp expansion"
echo "✓ Overlapping expanded intervals merge correctly"
echo "✓ Non-overlapping intervals remain separate"
echo ""

# Test 2: Verify chromosome-start clipping
cat > "$TEST_DIR/edge.bed" <<'EOF'
chr1	10	60
EOF

echo "Test 2: Chromosome-start clipping"
bedtools slop -i "$TEST_DIR/edge.bed" -g "$TEST_DIR/test.fai" -b 100 | \
  bedtools merge -i - > "$TEST_DIR/edge_padded.bed"

EDGE_LINE=$(cat "$TEST_DIR/edge_padded.bed")
# chr1:10-60 with 100bp padding should be chr1:0-160 (clipped at 0, not -90)
if [[ "$EDGE_LINE" != "chr1"$'\t'"0"$'\t'"160" ]]; then
  echo "FAIL: Chromosome-start clipping failed"
  echo "Expected: chr1 0 160"
  echo "Got: $EDGE_LINE"
  exit 1
fi

echo "✓ Reference bounds clipping works correctly"
echo ""

# Test 3: Verify mutect2.yaml uses --interval-padding 100
echo "Test 3: Mutect2 --interval-padding 100"
MUTECT2_PADDING_COUNT=$(/bin/grep -c "interval-padding 100" tasks/base_tasks/mutect2.yaml || echo 0)
if [[ "$MUTECT2_PADDING_COUNT" != "2" ]]; then
  echo "FAIL: mutect2.yaml should have exactly 2 occurrences of '--interval-padding 100'"
  echo "Found: $MUTECT2_PADDING_COUNT"
  exit 1
fi
echo "✓ Both Mutect2 commands use --interval-padding 100"
echo ""

# Test 4: Verify bcftools-filter uses padded BED in v2 workflows
echo "Test 4: Bcftools-filter uses validated-bed-padded-file"
for workflow in workflows/somatic/v2/somatic-dna-bam.yaml workflows/somatic/v2/somatic-dna-fastq.yaml; do
  # Check that post-mutect2 bcftools-filter uses validated-bed-padded-file
  if ! /bin/grep -q "validated-bed-padded-file" "$workflow"; then
    echo "FAIL: $workflow should use validated-bed-padded-file"
    exit 1
  fi
  echo "✓ $workflow bcftools-filter uses padded BED"
done
echo ""

# Test 5: Verify TMB still uses validated-bed-file (NOT padded)
echo "Test 5: TMB uses validated-bed-file (not padded)"
for workflow in workflows/somatic/v2/somatic-dna-bam.yaml workflows/somatic/v2/somatic-dna-fastq.yaml; do
  # Find collect-tmb-data task and verify it uses validated-bed-file (not padded)
  TMB_BED=$(/bin/grep -A20 "name: collect-tmb-data" "$workflow" | /bin/grep -o "validated-bed[^}]*" | head -1)
  if [[ "$TMB_BED" != "validated-bed-file" ]]; then
    echo "FAIL: $workflow TMB should use validated-bed-file, got: $TMB_BED"
    exit 1
  fi
  echo "✓ $workflow TMB correctly uses unpadded BED"
done
echo ""

# Test 6: Verify Mutect2 uses validated-bed-file (NOT padded)
echo "Test 6: Mutect2 uses validated-bed-file (not padded)"
for workflow in workflows/somatic/v2/somatic-dna-bam.yaml workflows/somatic/v2/somatic-dna-fastq.yaml; do
  # Find mutect2 task and verify it uses validated-bed-file (not padded)
  MUTECT2_BED=$(/bin/grep -A21 "name: mutect2$" "$workflow" | /bin/grep -A1 "regions-file" | tail -1 | /bin/grep -o "validated-bed[^}]*")
  if [[ "$MUTECT2_BED" != "validated-bed-file" ]]; then
    echo "FAIL: $workflow Mutect2 should use validated-bed-file, got: $MUTECT2_BED"
    exit 1
  fi
  echo "✓ $workflow Mutect2 correctly uses unpadded BED with --interval-padding"
done

echo ""
echo "=== All Tests Passed ==="
