#!/usr/bin/env python3

input_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/raw/mm9.ensGene.gtf"
output_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/genes-mm9.bed"

# Store chromosome, gene boundaries, and strand for each gene
genes = {}

with open(input_file, "r") as infile:
    for line in infile:
        # Skip GTF comment lines
        if line.startswith("#"):
            continue

        columns = line.strip().split("\t")

        # Use transcript entries to define the full gene span
        if columns[2] != "transcript":
            continue

        chrom = columns[0]
        start = int(columns[3])
        end = int(columns[4])
        strand = columns[6]
        attributes = columns[8]

        # Extract Ensembl gene ID from the attributes column
        gene_id = attributes.split('gene_id "')[1].split('"')[0]

        # Add a new gene or update its outermost transcript boundaries
        if gene_id not in genes:
            genes[gene_id] = [chrom, start, end, strand]
        else:
            if start < genes[gene_id][1]:
                genes[gene_id][1] = start
            if end > genes[gene_id][2]:
                genes[gene_id][2] = end

with open(output_file, "w") as outfile:
    for gene_id in genes:
        chrom = genes[gene_id][0]
        start = genes[gene_id][1] - 1
        end = genes[gene_id][2]
        strand = genes[gene_id][3]

        # Write gene-level BED6 output
        outfile.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n")