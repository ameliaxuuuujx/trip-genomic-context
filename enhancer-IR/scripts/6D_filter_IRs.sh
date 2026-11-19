#!/bin/bash
set -euo pipefail

BASE="/Users/cmdb/trip-genomic-context/enhancer-IR"
DATA="$BASE/data/processed"

IRS="$DATA/6C_intergenic-IRs-mm10.bed"
ENHANCERS="$DATA/6B_intergenic-enhancers-mm10.bed"

OUTPUT="$DATA/6D_filtered-IRs-mm10.bed"

bedtools window \
    -a "$IRS" \
    -b "$ENHANCERS" \
    -w 50000 \
    -v \
    > "$OUTPUT"

echo "Figure 6D enhancer filtering completed."

echo "Original intergenic IRs:"
wc -l < "$IRS"

echo "IRs retained (no enhancer within 50 kb):"
wc -l < "$OUTPUT"

echo "IRs excluded:"
awk 'END {print NR}' "$IRS" | {
    read total
    retained=$(wc -l < "$OUTPUT")
    echo $((total - retained))
}

echo "First 5 retained IRs:"
head -n 5 "$OUTPUT"
