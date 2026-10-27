#!/usr/bin/env python3

import csv

base = "/Users/cmdb/trip-genomic-context/enhancer-IR"

counts_file = base + "/data/processed/6B_gene-window-enhancer-counts-mm10.bed"
enhancer_file = base + "/data/processed/6B_intergenic-enhancers-mm10.bed"
background_file = base + "/data/processed/6B_intergenic-background-mm10.bed"
sizes_file = base + "/data/raw/mm10.chrom.sizes"

output_file = base + "/results/figure6B_summary.tsv"

# Read chromosome sizes
chrom_sizes = {}

with open(sizes_file) as f:
    for line in f:
        chrom, size = line.strip().split("\t")
        chrom_sizes[chrom] = int(size)

# Calculate intergenic genome length
intergenic_bp = 0

with open(background_file) as f:
    for line in f:
        columns = line.strip().split("\t")
        intergenic_bp += int(columns[2]) - int(columns[1])

# Count enhancers
enhancer_count = 0

with open(enhancer_file) as f:
    for line in f:
        if line.strip():
            enhancer_count += 1

# Expected enhancer density
background_density = enhancer_count / intergenic_bp

print("Intergenic enhancers:", enhancer_count)
print("Intergenic genome length:", intergenic_bp)
print("Background density:", background_density)

# Summarize enhancer counts by side and distance bin
results = {}

with open(counts_file) as f:
    for line in f:
        columns = line.strip().split("\t")

        chrom = columns[0]
        start = int(columns[1])
        end = int(columns[2])
        side = columns[4]
        distance_bin = int(columns[5])
        count = int(columns[6])

        if chrom not in chrom_sizes:
            continue

        # Correct window length at chromosome boundaries
        start = max(0, start)
        end = min(end, chrom_sizes[chrom])

        length = end - start

        if length <= 0:
            continue

        key = (side, distance_bin)

        if key not in results:
            results[key] = [0, 0, 0]

        results[key][0] += count
        results[key][1] += length
        results[key][2] += 1

# Calculate relative frequency
with open(output_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "side",
        "distance_bin",
        "distance_kb",
        "enhancer_count",
        "window_count",
        "observed_density",
        "relative_frequency"
    ])

    for side, distance_bin in sorted(results):
        count, total_bp, n_windows = results[(side, distance_bin)]

        observed_density = count / total_bp
        relative_frequency = observed_density / background_density
        distance_kb = (distance_bin - 0.5) * 10

        writer.writerow([
            side,
            distance_bin,
            distance_kb,
            count,
            n_windows,
            observed_density,
            relative_frequency
        ])

print("Output:", output_file)