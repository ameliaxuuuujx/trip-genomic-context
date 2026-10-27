#!/usr/bin/env python3

# Input file
gene_bed = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/genes-mm10-from-mm9.bed"

# Output file
output_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/6B_gene-windows-mm10.bed"

# Window settings
bin_size = 10000
max_distance = 200000

total_genes = 0
total_windows = 0

with open(gene_bed, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        columns = line.strip().split("\t")

        chrom = columns[0]
        start = int(columns[1])
        end = int(columns[2])
        gene_id = columns[3]
        strand = columns[5]

        total_genes += 1

        # Create windows on both sides of the gene
        for side in ["Start", "End"]:

            if strand == "+":
                boundary = start if side == "Start" else end
                direction = -1 if side == "Start" else 1
            elif strand == "-":
                boundary = end if side == "Start" else start
                direction = 1 if side == "Start" else -1
            else:
                continue

            for distance in range(0, max_distance, bin_size):
                if direction == 1:
                    window_start = boundary + distance
                    window_end = boundary + distance + bin_size
                else:
                    window_start = boundary - distance - bin_size
                    window_end = boundary - distance

                # Skip windows outside chromosome start
                if window_end <= 0:
                    continue

                window_start = max(0, window_start)

                outfile.write(
                    f"{chrom}\t{window_start}\t{window_end}\t"
                    f"{gene_id}\t{side}\t{distance // bin_size + 1}\n"
                )

                total_windows += 1

print("Genes processed:", total_genes)
print("Windows created:", total_windows)
print("Output:", output_file)