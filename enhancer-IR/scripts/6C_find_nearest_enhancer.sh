#!/bin/bash
set -euo pipefail

base="/Users/cmdb/trip-genomic-context/enhancer-IR"

irs="$base/data/processed/6C_IRs-away-from-TSS-mm10.bed"
enhancers="$base/data/processed/6B_intergenic-enhancers-mm10.bed"
output="$base/data/processed/6C_nearest-enhancer-mm10.bed"

# Sort IRs and enhancers by genomic coordinates
LC_ALL=C sort -k1,1 -k2,2n "$irs" > "$base/data/processed/6C_IRs-sorted-mm10.bed"
LC_ALL=C sort -k1,1 -k2,2n "$enhancers" > "$base/data/processed/6C_enhancers-sorted-mm10.bed"

# Find the nearest enhancer for each IR
bedtools closest \
    -a "$base/data/processed/6C_IRs-sorted-mm10.bed" \
    -b "$base/data/processed/6C_enhancers-sorted-mm10.bed" \
    -d \
    -t first \
    > "$output"

echo "IR-enhancer assignments:"
wc -l < "$output"

head -n 3 "$output"