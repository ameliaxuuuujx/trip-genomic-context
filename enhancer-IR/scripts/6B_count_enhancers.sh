#!/bin/bash

windows="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/6B_gene-windows-mm10.bed"
enhancers="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/6B_intergenic-enhancers-mm10.bed"
output="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/6B_gene-window-enhancer-counts-mm10.bed"

bedtools intersect \
    -a "$windows" \
    -b "$enhancers" \
    -c \
    > "$output"

echo "Windows with enhancer counts:"
wc -l < "$output"

head "$output"