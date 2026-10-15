#!/bin/bash

intergenic="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-intergenic-mm10.bed"
genes="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/genes-mm10-filtered-sorted.bed"
output="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-nearest-gene-filtered-mm10.bed"




# Find the nearest gene for each intergenic TRIP integration
bedtools closest -a "$intergenic" -b "$genes" -d > "$output"

# Check how many IRs are in the output
echo "Intergenic IRs with nearest gene:"
wc -l "$output"

# Preview the first few lines
head "$output"