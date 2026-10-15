# Goal 1: Reproduce Figure 6 from Akhtar et al. (2013)

## Step 1: Convert TRIP integration coordinates from mm9 to mm10

- Convert the original TRIP integration file (`mmc2.txt`) into BED format.
- Use UCSC LiftOver to convert TRIP integration coordinates from mm9 to mm10.
- Use R to merge the original TRIP count columns with the lifted mm10 coordinates by matching `pool_barcode`.
- A total of 17,857 mPGK reporter integration sites (IRs) were included in the processed dataset.

## Step 2: Reanalyze Figure 6A

### 2.1 Classify genes as expressed or non-expressed

- We used RNA-seq counts from 11 monoclonal mESC lines (`mmc5.txt`).
- A gene was classified as expressed if at least one sample had a count > 0; otherwise, it was classified as non-expressed.
- A total of 18,463 genes were classified as expressed and 3,906 as non-expressed.
- Because the RNA-seq count table uses gene IDs from the older mm9-era annotation, we lifted that gene annotation from mm9 to mm10 while retaining the original Ensembl gene IDs.
- 37,494 gene records were successfully converted; 187 failed conversion.

### 2.2 Identify intergenic TRIP IRs

- We identified intergenic IRs using the full gene annotation.
- A total of 8,472 intergenic IRs were identified.

### 2.3 Filter genes based on RNA-seq data

- We retained only genes with Ensembl IDs present in the RNA-seq count table.
- 22,241 of the 37,494 annotated genes were retained.

### 2.4 Find the nearest gene for each intergenic IR

- We recalculated the nearest gene using the filtered gene annotation.
- Therefore, the nearest gene refers to the nearest gene with available RNA-seq counts, which may not be the nearest gene in the full annotation.
- Nearest genes were assigned based on genomic distance using `bedtools closest`.
- Ties were retained, resulting in 8,475 IR-gene assignments for 8,472 intergenic IRs.

### 2.5 Assign gene expression status

- All 8,475 IR-gene assignments were successfully matched to RNA-seq expression status.
- No assignments were missing expression status.
- Among the 8,472 unique IRs, 7,857 were associated with expressed genes and 615 with non-expressed genes.

### 2.6 Assign gene side

- We classified IRs as being on the Start or End side of the nearest gene based on gene strand and genomic position.
- 4,778 IR-gene assignments were on the Start side, and 3,697 were on the End side.

### 2.7 Divide IRs into distance bins

- We divided IR-gene assignments into 10 equal-count distance bins separately for the Start and End sides.
- Expressed and non-expressed genes shared the same distance-bin boundaries on each side.
- All 8,475 IR-gene assignments were included in the analysis.

### 2.8 Calculate IR expression and reproduce Figure 6A

- We calculated reporter expression as the ratio of summed cDNA counts to summed gDNA counts across two technical replicates.
- For each TRIP pool, gene expression status, gene side, and distance bin, we calculated the median IR expression.
- We then calculated the mean and standard error of the mean (SEM) across the six mPGK pools.
- Of the 40 groups, 38 contained data from all six pools, while two contained data from five pools.
- We plotted IR expression as a function of distance from the nearest gene, separating the Start and End sides.

### 2.9 Limitations

- The original endogenous gene expression dataset used for Figure 6A could not be identified.
- We therefore used the RNA-seq count matrix from the 11 monoclonal mESC lines, which was provided for the Figure 6E analysis.
- Gene expression status was determined using raw RNA-seq counts rather than the original FPKM-based classification.
- Nearest genes were selected only from genes with available RNA-seq counts. This may change nearest-gene assignments and distances compared with the original analysis.
- The reporter expression normalization used in this reanalysis is simplified and may differ from the original TRIP normalization method.
- Our reanalysis included 7,857 expressed-associated IRs and 615 non-expressed-associated IRs, compared with 5,435 and 2,746, respectively, in the published Figure 6A.
- Therefore, our current analysis reproduces the general relationship between gene proximity and reporter expression but does not fully reproduce the original Figure 6A.

### 2.10 Output files

- `data/processed/trip-distance-bins-filtered-mm10.tsv`: Processed IR-gene assignments with expression status, gene side, and distance bins.
- `results/figure6A_summary.tsv`: Mean IR expression, SEM, and distance information for each group.
- `results/figure6A_280kb.png`: Figure 6A reanalysis plot.

## Step 3: Reanalyze Figure 6B

To be completed.