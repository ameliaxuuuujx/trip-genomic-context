#!/usr/bin/env python3

import csv
import subprocess
from pathlib import Path

BASE = Path("/Users/cmdb/trip-genomic-context/enhancer-IR")
DATA = BASE / "data/processed"

OBSERVED = DATA / "6B_enhancer-nearest-boundary-mm10.bed"
RANDOM = DATA / "6B_random-nearest-boundary-mm10.bed"
GENES = DATA / "genes-mm10-from-mm9.bed"

INTERGENIC_RANDOM = DATA / "6B_intergenic-random-nearest-boundary-mm10.bed"
OUTPUT = DATA / "6B_relative-enhancer-frequency-mm10.tsv"

BIN_SIZE = 10000
MAX_DISTANCE = 200000


# Step 1: Keep only intergenic random positions
with open(INTERGENIC_RANDOM, "w") as out:
    subprocess.run(
        [
            "bedtools", "intersect",
            "-a", str(RANDOM),
            "-b", str(GENES),
            "-v"
        ],
        stdout=out,
        check=True
    )


# Step 2: Read nearest-boundary distances
def read_distances(filename):
    distances = {"Start": [], "End": []}

    with open(filename) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")

            side = cols[-2]
            distance = int(cols[-1])

            if side in distances and distance >= 0:
                distances[side].append(distance)

    return distances


observed = read_distances(OBSERVED)
random_all = read_distances(RANDOM)
random_intergenic = read_distances(INTERGENIC_RANDOM)


# Step 3: Count positions in distance bins
def count_bins(values):
    counts = [0] * (MAX_DISTANCE // BIN_SIZE)

    for distance in values:
        if 0 <= distance < MAX_DISTANCE:
            bin_index = distance // BIN_SIZE
            counts[bin_index] += 1

    return counts


# Step 4: Calculate observed/expected frequency
results = []

backgrounds = {
    "genome_wide": random_all,
    "intergenic_only": random_intergenic
}

for background_name, background in backgrounds.items():
    for side in ["Start", "End"]:

        obs = observed[side]
        exp = background[side]

        obs_counts = count_bins(obs)
        exp_counts = count_bins(exp)

        for i in range(len(obs_counts)):
            start = i * BIN_SIZE
            end = start + BIN_SIZE

            observed_frequency = obs_counts[i] / len(obs)
            expected_frequency = exp_counts[i] / len(exp)

            if expected_frequency > 0:
                relative_frequency = (
                    observed_frequency / expected_frequency
                )
            else:
                relative_frequency = float("nan")

            results.append({
                "background": background_name,
                "side": side,
                "distance_start_bp": start,
                "distance_end_bp": end,
                "distance_midpoint_kb": (start + end) / 2000,
                "observed_count": obs_counts[i],
                "expected_count": exp_counts[i],
                "observed_frequency": observed_frequency,
                "expected_frequency": expected_frequency,
                "relative_frequency": relative_frequency
            })


# Step 5: Save results
with open(OUTPUT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(results)


# Step 6: Print summary
print("Figure 6B relative frequency calculation completed.")

for side in ["Start", "End"]:
    print(f"Observed {side}: {len(observed[side])}")
    print(f"Random genome-wide {side}: {len(random_all[side])}")
    print(f"Random intergenic {side}: {len(random_intergenic[side])}")

print("Total random positions:", sum(map(len, random_all.values())))
print(
    "Intergenic random positions:",
    sum(map(len, random_intergenic.values()))
)

print("Distance bins per side:", MAX_DISTANCE // BIN_SIZE)
print("Output:", OUTPUT)
