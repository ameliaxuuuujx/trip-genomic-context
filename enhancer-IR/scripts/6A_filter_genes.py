#!/usr/bin/env python3

# Input files
counts_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/raw/mmc5.txt"
gene_bed = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/genes-mm10-from-mm9.bed"

# Output file
output_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/genes-mm10-filtered.bed"

# Get gene IDs from RNA-seq counts
gene_ids = set()

with open(counts_file, "r") as f:
    next(f)  # Skip header

    for line in f:
        columns = line.strip().split("\t")
        gene_ids.add(columns[0])

print("Genes in RNA-seq:", len(gene_ids))


# Keep genes with RNA-seq counts
total = 0
matched = 0

with open(gene_bed, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        columns = line.strip().split("\t")
        gene_id = columns[3]

        total += 1

        if gene_id in gene_ids:
            outfile.write(line)
            matched += 1

print("Genes in BED:", total)
print("Genes kept:", matched)
print("Genes removed:", total - matched)