#!/bin/bash

base="/Users/cmdb/trip-genomic-context/enhancer-IR"

genes="$base/data/processed/genes-mm10-from-mm9.bed"
sizes="$base/data/raw/mm10.chrom.sizes"
output="$base/data/processed/6B_intergenic-background-mm10.bed"

# Download mm10 chromosome sizes
curl -fL \
  "https://hgdownload.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes" \
  -o "$sizes" || exit 1

# Check downloaded file
if [ ! -s "$sizes" ]; then
    echo "Error: chromosome sizes file is empty."
    exit 1
fi

# Sort and merge genes, then find intergenic regions
bedtools sort -i "$genes" -g "$sizes" |
    bedtools merge -i - |
    bedtools complement -i - -g "$sizes" > "$output"

echo "Intergenic background:"
wc -l < "$output"