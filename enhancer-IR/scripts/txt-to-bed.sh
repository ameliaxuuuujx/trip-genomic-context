#!/bin/bash


input="/Users/cmdb/trip-genomic-context/enhancer-IR/data/raw/mmc2.txt"
output="/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-mm9-short.bed"

awk 'BEGIN {OFS="\t"}
     /^#/ {next}
     $1=="TRIP_pool" {next}
     {
       print $3, $4-1, $4, $1"_"$2, 0, $5
     }' "$input" > "$output"

echo "Created $output"