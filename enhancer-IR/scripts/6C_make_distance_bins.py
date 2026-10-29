#!/usr/bin/env python3

import csv

base = "/Users/cmdb/trip-genomic-context/enhancer-IR"

input_file = base + "/data/processed/6C_nearest-enhancer-mm10.bed"
output_file = base + "/data/processed/6C_distance-bins-mm10.tsv"

# Read IR-enhancer assignments and exclude invalid distances
rows = []
excluded = 0

with open(input_file) as f:
    reader = csv.reader(f, delimiter="\t")

    for cols in reader:
        if not cols:
            continue

        if len(cols) != 16:
            raise ValueError(f"Expected 16 columns, found {len(cols)}")

        if int(cols[15]) == -1:
            excluded += 1
            continue

        rows.append(cols)

# Sort IRs by distance to the nearest enhancer
rows.sort(key=lambda row: int(row[15]))

# Assign 10 approximately equal-sized distance bins
n = len(rows)

if n == 0:
    raise ValueError("No valid IRs found")

for i, row in enumerate(rows):
    distance_bin = i * 10 // n + 1
    row.append(str(distance_bin))

# Write the output with column names
columns = [
    "chrom",
    "start",
    "end",
    "pool_barcode",
    "score",
    "strand",
    "normalization_counts_1",
    "normalization_counts_2",
    "expression_counts_1",
    "expression_counts_2",
    "enhancer_chrom",
    "enhancer_start",
    "enhancer_end",
    "enhancer_id",
    "enhancer_category",
    "distance",
    "distance_bin"
]

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(columns)
    writer.writerows(rows)

# Print summary
print("Excluded IRs without matching enhancers:", excluded)
print("Total valid IRs:", n)
print("Distance bins:", 10)

for bin_number in range(1, 11):
    bin_rows = [
        row for row in rows
        if int(row[16]) == bin_number
    ]

    distances = [int(row[15]) for row in bin_rows]

    print(
        f"Bin {bin_number}: "
        f"{len(bin_rows)} IRs, "
        f"distance {min(distances)}–{max(distances)} bp"
    )

print("Output:", output_file)
