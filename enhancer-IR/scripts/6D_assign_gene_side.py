#!/usr/bin/env python3

BASE = "/Users/cmdb/trip-genomic-context/enhancer-IR"

INPUT = BASE + "/data/processed/6D_nearest-gene-mm10.bed"
OUTPUT = BASE + "/data/processed/6D_gene-side-mm10.tsv"

start_count = 0
end_count = 0

with open(INPUT) as f, open(OUTPUT, "w") as out:

    for line in f:
        cols = line.rstrip("\n").split("\t")

        if len(cols) != 17:
            raise ValueError(f"Expected 17 columns, got {len(cols)}")

        ir_start = int(cols[1])
        ir_end = int(cols[2])

        gene_start = int(cols[11])
        gene_end = int(cols[12])
        gene_strand = cols[15]

        if ir_end <= gene_start:
            side = "Start" if gene_strand == "+" else "End"

        elif ir_start >= gene_end:
            side = "End" if gene_strand == "+" else "Start"

        else:
            raise ValueError("IR overlaps nearest gene")

        if side == "Start":
            start_count += 1
        else:
            end_count += 1

        out.write("\t".join(cols + [side]) + "\n")

print("Figure 6D gene-side assignment completed.")
print("Start side:", start_count)
print("End side:", end_count)
print("Total assignments:", start_count + end_count)
print("Output:", OUTPUT)
