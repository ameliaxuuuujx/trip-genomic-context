#!/usr/bin/env python3

# Input file
input_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/raw/ENCFF110IIC.bed"

# Output file
output_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/6B_enhancers-mm10.bed"

# Filter enhancer-like elements
total = 0
matched = 0

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        columns = line.strip().split("\t")
        total += 1

        category = columns[9]

        if category == "dELS" or category == "pELS":
            outfile.write("\t".join([
                columns[0],  # Chromosome
                columns[1],  # Start
                columns[2],  # End
                columns[3],  # Element ID
                category     # Enhancer category
            ]) + "\n")

            matched += 1

print("Total cCREs:", total)
print("Enhancer-like elements:", matched)
print("Other elements removed:", total - matched)