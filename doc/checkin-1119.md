# Project Check-in — November 19, 2026

## 1. Addressing Previous Feedback

We revised Figure 6B using one million random genomic coordinates as the background instead of genome-wide enhancer density. We also applied quantile-based binning to Figure 6C, as suggested. Our scripts remain modular, and we plan to add a master script to automate the pipeline.

## 2. New Progress

**Figure 6B:** We completed the revised enhancer enrichment analysis using both genome-wide and intergenic-only random backgrounds. Both showed enrichment near gene boundaries, resolving the issue of all bins being above background.

**Figure 6C:** We completed the enhancer proximity analysis using 3,058 IRs and 10 quantile bins. We also added a 0–50 kb analysis to examine short-range effects.

**Figure 6D:** We completed the gene proximity analysis using 5,636 IRs after excluding those with nearby enhancers. Reporter expression increased near gene boundaries, consistent with the original paper. We also added a 0–50 kb analysis.

**Figure 6E:** [To be completed by Samantha.]

## 3. Project Organization

We maintain separate `enhancer-IR/` and `rnaseq/` directories, with modular scripts, detailed READMEs, and summary tables for reproducibility. Large intermediate files are excluded from GitHub.

## 4. Challenges and Questions

- **Background selection:** Which background is more appropriate for Figure 6B: genome-wide or intergenic-only random coordinates?
- **Sample sizes:** Figures 6C and 6D include more IRs than the original study because of differences in enhancer annotations and filtering. How should we best address these differences?
- **Pipeline:** We plan to integrate our scripts into a master pipeline.
- **Figure 6E:** [To be completed by Samantha.]