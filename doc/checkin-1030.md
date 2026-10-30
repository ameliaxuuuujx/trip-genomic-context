# Checkin-1030

## 1. Addressing Previous Feedback

We received feedback on our project proposal regarding code modularity, the use of existing bioinformatics tools, and computational resource requirements. We have addressed these suggestions in the following ways:

- **Code modularity**: We organized our analysis into separate scripts for different processing steps, allowing us to modify individual steps without rerunning the entire pipeline and making the analysis easier to reproduce.

- **Use of existing tools**: Rather than implementing genomic interval operations from scratch, we used established tools such as bedtools for identifying intergenic regions, finding nearest genomic features, and counting overlaps. We also used Python for data processing and R with tidyverse for statistical summaries and visualization.

- **Computational resources**: We are working with publicly available processed TRIP data rather than raw sequencing reads, and our analyses so far have been successfully performed on personal laptops without requiring HPC resources. For our planned Micro-C extension, we have also identified an appropriate software platform for analyzing chromatin contact data, so we do not currently anticipate needing HPC resources for this part of the project either.

## 2. Progress Since Proposal

- **Data preparation**: We obtained the published TRIP datasets and converted the integration coordinates from mm9 to mm10 using UCSC LiftOver. We successfully processed 17,857 mPGK reporter integration sites (IRs) for downstream analysis.

- **Figure 6A – Gene proximity and reporter expression**: We identified 8,472 intergenic IRs and analyzed their expression relative to neighboring endogenous genes. We classified genes as expressed or non-expressed using available RNA-seq data and calculated median reporter expression across distance bins. Our reanalysis reproduced the general relationship between gene proximity and reporter expression, although the proportions of IRs associated with expressed and non-expressed genes differed from the original study.

- **Figure 6B – Enhancer distribution around genes**: We used ENCODE cCRE annotations to identify 7,871 intergenic enhancer-like elements. We counted their occurrence in 10 kb windows extending up to 200 kb from gene boundaries and calculated relative enhancer frequencies using genome-wide intergenic enhancer density as the background.

- **Figure 6C – Enhancer proximity and reporter expression (in progress)**: We have begun analyzing how reporter expression varies with distance from the nearest enhancer, while excluding IRs located within 50 kb of endogenous gene transcription start sites. We are currently exploring different distance-binning strategies to better characterize short-range enhancer effects.

- **Figure 6E: Effect of IRs on endogenous gene expression**: We ran DESeq2 on the RNA-seq clones and applied the holms method to get adjusted values of significance to account for multiple testing. WE used bedtools window to extract the genes within 100kb of an IR, and then we used R studio to calculate upstream and downstream distance. One challenge is the plot does not match the figure, which may be due to issues with statistical testing.

## 3. Project Organization

We organized our GitHub repository into separate directories for different analyses. Figures 6A–6D are grouped under enhancer-IR/, while Figure 6E is organized under rnaseq/. Each analysis directory contains data, scripts, and results, with a README documenting the workflow and outputs.

Our code is modular, with separate scripts for data processing, genomic feature analysis, and visualization. Large intermediate files are excluded from GitHub using .gitignore but can be regenerated using the provided scripts.

We divided the work between team members, with Amelia focusing on Figures 6A–6D and Samantha working on Figure 6E. We use GitHub to share code, track progress, and integrate our results.

## 4. Challenges and Questions

- **Figrue 6A**: We could not identify the exact endogenous gene expression dataset used in the original Figure 6A. Instead, we used RNA-seq data from 11 monoclonal mESC lines provided for Figure 6E. This resulted in different proportions of expressed and non-expressed genes compared with the published results. We would appreciate advice on whether this is an appropriate alternative or whether another dataset would be more suitable.

- **Figure 6B**: In our reanalysis, the relative enhancer frequency remained above the background level (y = 1) across all distance bins, even at larger distances from genes. We used genome-wide intergenic enhancer density as the background, but this may not accurately represent the expected enhancer frequency around genes. We would appreciate advice on whether a random genomic background would be more appropriate.

- **Figure 6C (ongoing)**: Enhancer distances are unevenly distributed, making it challenging to visualize short-range effects. We are exploring equal-count bins rather than fixed-width bins and would appreciate advice on the most appropriate binning strategy.

- **Figure 6E (ongoing)**: We have a plot, but it does not reflect the trend shown in the figure. While the figure shows a clear relationship between gene proximity to IR and differential expression, this was not clear in the reproduced figure. Notably, two observations are that there are much fewer data points as well as significant genes plotted. To address this, the next step is to look where these points were lost, which we expect could be due to the GTF file or due to the multiple hypothesis testing step. The data could have been lost when adding coordinates from the GTF file, we lost 1,614 genes (7.14%), and many genes lost significance after accounting for multiple hypothesis testing. In addition to looking deeper into these potential issues, we will look at visualize several locus specific examples to gain more confidence in the genome-wide observations.
