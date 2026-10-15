#!/bin/bash

trip="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-mm10-long.bed"
genes="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/genes-mm10-from-mm9.bed"

trip_sorted="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-mm10-long-sorted.bed"
genes_sorted="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/genes-mm10-sorted.bed"
output="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-intergenic-mm10.bed"

# Sort TRIP integrations and genes by chromosome and start position
LC_ALL=C sort -k1,1 -k2,2n "$trip" > "$trip_sorted"
LC_ALL=C sort -k1,1 -k2,2n "$genes" > "$genes_sorted"

# Keep only TRIP integrations that do not overlap any gene
bedtools intersect -a "$trip_sorted" -b "$genes_sorted" -v > "$output"

# Show how many IRs were present before and after filtering
echo "All TRIP IRs:"
wc -l "$trip_sorted"
### All TRIP IRs: 17857

echo "Intergenic TRIP IRs:"
wc -l "$output"
### Intergenic TRIP IRs: 8210