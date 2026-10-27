#!/bin/bash

enhancers="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/6B_enhancers-mm10.bed"
genes="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/genes-mm10-from-mm9.bed"
output="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/6B_intergenic-enhancers-mm10.bed"

bedtools intersect \
    -a "$enhancers" \
    -b "$genes" \
    -v \
    > "$output"

echo "Intergenic enhancer-like elements:"
wc -l < "$output"