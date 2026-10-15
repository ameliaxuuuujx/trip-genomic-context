# Goal 1

## Step 1: Convert TRIP integration coordinates from mm9 to mm10
- Convert the original TRIP integration file (`mmc2.txt`) into BED format.
- Use UCSC LiftOver to convert TRIP integration coordinates from mm9 to mm10.
- Use R to merge the original TRIP count columns with the lifted mm10 coordinates by matching `pool_barcode`.

## Step 2: Reanalyze Figure 6A

- **Classify genes as expressed or non-expressed**
    - We used RNA-seq counts from 11 monoclonal mESC lines.
    - A gene was classified as expressed if at least one sample had a count > 0; otherwise, it was classified as non-expressed.YUIO
    - Because the RNA-seq count table uses gene IDs from the older mm9-era annotation, we lifted that gene annotation from mm9 to mm10 while retaining the original gene IDs.
    - 37,494 records were successfully converted; 187 failed conversion.

- **Identify intergenic TRIP IRs**
    - We identified 8,472 intergenic TRIP IRs.

- **Find the nearest gene for each intergenic IR**
    - Nearest genes were assigned based on genomic distance.
    - Ties were retained; 6 IRs had equally close nearest genes.

- **Assign gene expression status**
    - 5,437 IR-gene assignments matched genes in the RNA-seq count table. not expressed + expressed
    - 3,041 assignments could not be matched and were treated as missing for expression-status analysis.  ABCDE