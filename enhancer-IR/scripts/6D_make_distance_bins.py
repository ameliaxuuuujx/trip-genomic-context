#!/usr/bin/env python3

import csv

BASE = "/Users/cmdb/trip-genomic-context/enhancer-IR"

INPUT = BASE + "/data/processed/6D_gene-side-mm10.tsv"
OUTPUT = BASE + "/data/processed/6D_distance-bins-mm10.tsv"

columns = [
    "ir_chrom", "ir_start", "ir_end", "pool_barcode",
    "ir_score", "ir_strand",
    "normalization_counts_1", "normalization_counts_2",
    "expression_counts_1", "expression_counts_2",
    "gene_chrom", "gene_start", "gene_end",
    "gene_id", "gene_score", "gene_strand",
    "distance", "side"
]

# Read IR-gene assignments
data = {"Start": [], "End": []}

with open(INPUT) as f:
    reader = csv.reader(f, delimiter="\t")

    for row in reader:
        if len(row) != 18:
            raise ValueError(f"Expected 18 columns, got {len(row)}")

        side = row[17]
        distance = int(row[16])

        if side not in data:
            raise ValueError(f"Unknown gene side: {side}")

        if distance < 0:
            raise ValueError(f"Invalid distance: {distance}")

        data[side].append(row)

# Assign 10 equal-count bins separately for each side
results = []

for side in ["Start", "End"]:
    rows = sorted(data[side], key=lambda row: int(row[16]))
    n = len(rows)

    for i, row in enumerate(rows):
        distance_bin = (i * 10 // n) + 1
        results.append(row + [distance_bin])

    print(f"\n{side} side: {n} IRs")

    for bin_number in range(1, 11):
        bin_rows = [
            row for row in results
            if row[17] == side and row[18] == bin_number
        ]

        distances = [int(row[16]) for row in bin_rows]

        print(
            f"Bin {bin_number}: {len(bin_rows)} IRs, "
            f"distance {min(distances)}–{max(distances)} bp"
        )

# Save results with header
with open(OUTPUT, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(columns + ["distance_bin"])
    writer.writerows(results)

print("\nFigure 6D distance binning completed.")
print("Total IRs:", len(results))
print("Output:", OUTPUT)
