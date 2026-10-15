#!/usr/bin/env python3

status_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/gene-expression-status.tsv"
nearest_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-nearest-gene-filtered-mm10.bed"
output_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-nearest-gene-status-filtered-mm10.bed"

# Store gene expression status in a dictionary
gene_status = {}

with open(status_file, "r") as infile:
    # Skip header
    infile.readline()

    for line in infile:
        columns = line.strip().split("\t")

        gene_id = columns[0]
        status = columns[1]

        gene_status[gene_id] = status

matched = 0
missing = 0

with open(nearest_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        columns = line.strip().split("\t")

        # Gene ID is column 14 in this file
        gene_id = columns[13]

        if gene_id in gene_status:
            status = gene_status[gene_id]
            matched += 1
        else:
            status = "NA"
            missing += 1

        outfile.write(line.strip() + "\t" + status + "\n")

print("Matched genes:", matched)
print("Missing genes:", missing)
print("Output:", output_file)