#!/usr/bin/env python3

input_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-gene-side-filtered-mm10.bed"
output_file = "/Users/cmdb/trip-genomic-context/enhancer-IR/data/processed/trip-distance-bins-filtered-mm10.tsv"

start_rows = []
end_rows = []

with open(input_file, "r") as infile:
    for line in infile:
        columns = line.strip().split("\t")

        distance = int(columns[16])
        status = columns[17]
        side = columns[18]

        # Exclude genes without RNA-seq expression status
        if status == "NA":
            continue

        if side == "Start":
            start_rows.append((distance, columns))

        elif side == "End":
            end_rows.append((distance, columns))


# Sort IRs from nearest to farthest
start_rows.sort(key=lambda x: x[0])
end_rows.sort(key=lambda x: x[0])


def assign_bins(rows):
    total = len(rows)
    output = []

    for i in range(total):
        distance = rows[i][0]
        columns = rows[i][1]

        # Divide rows into 10 approximately equal-sized bins
        bin_number = int(i * 10 / total) + 1

        output.append((columns, bin_number))

    return output


start_binned = assign_bins(start_rows)
end_binned = assign_bins(end_rows)


with open(output_file, "w") as outfile:
    for columns, bin_number in start_binned:
        outfile.write("\t".join(columns) + "\t" + str(bin_number) + "\n")

    for columns, bin_number in end_binned:
        outfile.write("\t".join(columns) + "\t" + str(bin_number) + "\n")


print("Start-side IR-gene rows:", len(start_rows))
print("End-side IR-gene rows:", len(end_rows))
print("Total rows used:", len(start_rows) + len(end_rows))
print("Output:", output_file)