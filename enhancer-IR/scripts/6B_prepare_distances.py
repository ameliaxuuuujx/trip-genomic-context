#!/usr/bin/env python3

base = "/Users/cmdb/trip-genomic-context/enhancer-IR"

gene_file = base + "/data/processed/genes-mm10-from-mm9.bed"
enhancer_file = base + "/data/processed/6B_intergenic-enhancers-mm10.bed"

gene_output = base + "/data/processed/6B_gene-boundaries-mm10.bed"
enhancer_output = base + "/data/processed/6B_enhancer-midpoints-mm10.bed"

# Extract gene Start and End positions
gene_count = 0

with open(gene_file) as f, open(gene_output, "w") as out:
    for line in f:
        cols = line.strip().split("\t")

        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        gene_id = cols[3]
        strand = cols[5]

        if strand == "+":
            gene_start = start
            gene_end = end - 1
        elif strand == "-":
            gene_start = end - 1
            gene_end = start
        else:
            continue

        out.write(
            f"{chrom}\t{gene_start}\t{gene_start + 1}\t{gene_id}\tStart\n"
        )
        out.write(
            f"{chrom}\t{gene_end}\t{gene_end + 1}\t{gene_id}\tEnd\n"
        )

        gene_count += 1

# Extract enhancer midpoints
enhancer_count = 0

with open(enhancer_file) as f, open(enhancer_output, "w") as out:
    for line in f:
        cols = line.strip().split("\t")

        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        enhancer_id = cols[3]

        midpoint = (start + end) // 2

        out.write(
            f"{chrom}\t{midpoint}\t{midpoint + 1}\t{enhancer_id}\n"
        )

        enhancer_count += 1

print("Genes processed:", gene_count)
print("Gene boundaries:", gene_count * 2)
print("Enhancer midpoints:", enhancer_count)
print("Output:", gene_output)
print("Output:", enhancer_output)
