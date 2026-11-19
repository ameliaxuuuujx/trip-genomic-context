#!/bin/bash
set -euo pipefail

BASE="/Users/cmdb/trip-genomic-context/enhancer-IR"

INPUT="$BASE/data/processed"
OUTPUT="$BASE/data/processed"

bedtools sort \
    -i "$INPUT/6B_gene-boundaries-mm10.bed" \
    > "$OUTPUT/6B_gene-boundaries-sorted-mm10.bed"

bedtools sort \
    -i "$INPUT/6B_enhancer-midpoints-mm10.bed" \
    > "$OUTPUT/6B_enhancer-midpoints-sorted-mm10.bed"

bedtools closest \
    -a "$OUTPUT/6B_enhancer-midpoints-sorted-mm10.bed" \
    -b "$OUTPUT/6B_gene-boundaries-sorted-mm10.bed" \
    -d -t first \
    > "$OUTPUT/6B_enhancer-nearest-boundary-mm10.bed"

echo "Enhancer-to-nearest-boundary calculation completed."
echo "Number of enhancer assignments:"
wc -l < "$OUTPUT/6B_enhancer-nearest-boundary-mm10.bed"

echo "First 5 rows:"
head -n 5 "$OUTPUT/6B_enhancer-nearest-boundary-mm10.bed"
