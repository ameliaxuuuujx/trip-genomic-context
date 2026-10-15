#!/usr/bin/env python3

input_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/raw/mmc5.txt"
output_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/gene-expression-status.tsv"

expressed_count = 0
nonexpressed_count = 0

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    # Skip the header line
    infile.readline()

    # Write the output header
    outfile.write("gene_id\tstatus\n")

    for line in infile:
        columns = line.strip().split("\t")

        gene_id = columns[0]

        # Columns 2-12 are RNA-seq counts from the 11 clones
        counts = []

        for value in columns[1:]:
            counts.append(int(value))

        # A gene is expressed if at least one clone has a non-zero count
        if max(counts) > 0:
            status = "expressed"
            expressed_count += 1
        else:
            status = "non-expressed"
            nonexpressed_count += 1

        outfile.write(f"{gene_id}\t{status}\n")

print("Expressed genes:", expressed_count)
print("Non-expressed genes:", nonexpressed_count)
print("Output:", output_file)