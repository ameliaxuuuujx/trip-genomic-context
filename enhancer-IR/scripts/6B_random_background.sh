#!/bin/bash
set -euo pipefail

BASE="/Users/cmdb/trip-genomic-context/enhancer-IR"
DATA="$BASE/data/processed"

GENOME="$BASE/data/raw/mm10.chrom.sizes"
BOUNDARIES="$DATA/6B_gene-boundaries-sorted-mm10.bed"

RANDOM_FILE="$DATA/6B_random-positions-mm10.bed"
SORTED="$DATA/6B_random-positions-sorted-mm10.bed"
OUTPUT="$DATA/6B_random-nearest-boundary-mm10.bed"

bedtools random \
    -l 1 \
    -n 1000000 \
    -g "$GENOME" \
    -seed 42 \
    > "$RANDOM_FILE"

bedtools sort \
    -i "$RANDOM_FILE" \
    > "$SORTED"

bedtools closest \
    -a "$SORTED" \
    -b "$BOUNDARIES" \
    -d -t first \
    > "$OUTPUT"

echo "Random background calculation completed."

echo "Number of random positions:"
wc -l < "$OUTPUT"

echo "First 5 rows:"
head -n 5 "$OUTPUT"
