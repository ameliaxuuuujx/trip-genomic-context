#!/usr/bin/env bash

#Use bedtools window to extract all genes within 100kb of each IR
#Make a loop to run this for every clone_rnaseq to every clone_IRs.

cd /Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles

#bedtools window -a {clone}_IR.bed -b {clone}genes.bed -w 100000 

for clone in 16 17 18 19 20 21 22 23 24 25 26
do
  bedtools window \
    -a /Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR_${clone}.bed \
    -b /Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/C${clone}genes.bed \
    -w 100000 \
    > /Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR${clone}_window.bed

  echo "Done IR${clone}"
done
