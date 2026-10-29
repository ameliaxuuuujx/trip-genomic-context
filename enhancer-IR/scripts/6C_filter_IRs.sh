#!/bin/bash
set -euo pipefail

base="/Users/cmdb/trip-genomic-context/enhancer-IR"

irs="$base/data/processed/6C_intergenic-IRs-mm10.bed"
tss="$base/data/processed/6C_gene-TSS-mm10.bed"

output="$base/data/processed/6C_IRs-away-from-TSS-mm10.bed"

# Exclude IRs within 50 kb of any endogenous gene TSS
bedtools window \
    -a "$irs" \
    -b "$tss" \
    -w 50000 \
    -v > "$output"

echo "Intergenic IRs before filtering:"
wc -l < "$irs"

echo "IRs remaining after excluding nearby TSSs:"
wc -l < "$output"