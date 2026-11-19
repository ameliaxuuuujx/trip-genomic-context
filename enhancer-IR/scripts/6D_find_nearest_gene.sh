#!/bin/bash
set -euo pipefail

BASE="/Users/cmdb/trip-genomic-context/enhancer-IR"
DATA="$BASE/data/processed"

IRS="$DATA/6D_filtered-IRs-mm10.bed"
GENES="$DATA/genes-mm10-from-mm9.bed"

SORTED_IRS="$DATA/6D_filtered-IRs-sorted-mm10.bed"
SORTED_GENES="$DATA/6D_genes-sorted-mm10.bed"

OUTPUT="$DATA/6D_nearest-gene-mm10.bed"

bedtools sort \
    -i "$IRS" \
    > "$SORTED_IRS"

bedtools sort \
    -i "$GENES" \
    > "$SORTED_GENES"

bedtools closest \
    -a "$SORTED_IRS" \
    -b "$SORTED_GENES" \
    -d -t first \
    > "$OUTPUT"

echo "Figure 6D nearest-gene calculation completed."

echo "Number of IR-gene assignments:"
wc -l < "$OUTPUT"

echo "First 5 rows:"
head -n 5 "$OUTPUT"
