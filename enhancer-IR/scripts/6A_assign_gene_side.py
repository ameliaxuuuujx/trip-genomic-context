#!/usr/bin/env python3

input_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-nearest-gene-status-mm10.bed"
output_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-gene-side-mm10.bed"

start_side = 0
end_side = 0
missing_status = 0

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        columns = line.strip().split("\t")

        ir_start = int(columns[1])
        ir_end = int(columns[2])

        gene_start = int(columns[11])
        gene_end = int(columns[12])
        gene_strand = columns[15]

        status = columns[17]

        # Keep track of rows without expression status
        if status == "NA":
            missing_status += 1

        # IR is on the left side of the gene
        if ir_end <= gene_start:
            if gene_strand == "+":
                side = "Start"
            else:
                side = "End"

        # IR is on the right side of the gene
        elif ir_start >= gene_end:
            if gene_strand == "+":
                side = "End"
            else:
                side = "Start"

        else:
            side = "overlap"

        if side == "Start":
            start_side += 1
        elif side == "End":
            end_side += 1

        outfile.write(line.strip() + "\t" + side + "\n")

print("Start-side IRs:", start_side)
print("End-side IRs:", end_side)
print("IRs with missing expression status:", missing_status)
print("Output:", output_file)