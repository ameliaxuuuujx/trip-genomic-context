#!/usr/bin/env python3

import csv

base = "/Users/cmdb/trip-genomic-context/enhancer-IR"

ir_file = base + "/data/processed/trip-distance-bins-filtered-mm10.tsv"
gene_file = base + "/data/processed/genes-mm10-from-mm9.bed"

ir_output = base + "/data/processed/6C_intergenic-IRs-mm10.bed"
tss_output = base + "/data/processed/6C_gene-TSS-mm10.bed"

# Extract unique intergenic IRs from Figure 6A data
seen = set()

with open(ir_file) as f, open(ir_output, "w") as out:
    # Figure 6A TSV has no header
    columns = [
        "chrom", "start", "end", "pool_barcode", "score", "strand",
        "normalization_counts_1", "normalization_counts_2",
        "expression_counts_1", "expression_counts_2",
        "gene_chrom", "gene_start", "gene_end", "gene_id",
        "gene_score", "gene_strand", "distance", "gene_status",
        "side", "distance_bin"
    ]

    reader = csv.DictReader(f, fieldnames=columns, delimiter="\t")

    for row in reader:
        barcode = row["pool_barcode"]

        if barcode in seen:
            continue

        seen.add(barcode)

        out.write("\t".join([
            row["chrom"],
            row["start"],
            row["end"],
            barcode,
            row["score"],
            row["strand"],
            row["normalization_counts_1"],
            row["normalization_counts_2"],
            row["expression_counts_1"],
            row["expression_counts_2"]
        ]) + "\n")

print("Unique intergenic IRs:", len(seen))

# Extract transcription start sites from gene annotations
tss_count = 0

with open(gene_file) as f, open(tss_output, "w") as out:
    for line in f:
        cols = line.strip().split("\t")

        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        gene_id = cols[3]
        strand = cols[5]

        # TSS depends on gene orientation
        if strand == "+":
            tss = start
        elif strand == "-":
            tss = end - 1
        else:
            continue

        out.write(
            f"{chrom}\t{tss}\t{tss + 1}\t{gene_id}\t0\t{strand}\n"
        )

        tss_count += 1

print("Gene TSS records:", tss_count)