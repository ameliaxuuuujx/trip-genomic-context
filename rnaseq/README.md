Differential analysis of 11 TRIP Clones from the raw counts.

The summarized output (will be) in `rnaseq/results/results.rmd`

A TRIP clone are single cells selected from the transgene pool in order to test the impact of IRs on transcription. Each of the 11 clones has their own unique integrations in random locations throughout the genome. To run differnential analysis of the rna counts, the authors contrast each clone to the remaining 10. Here are the steps:

- 261009: Created metadata tables containing the clone name and their group (exp= the clone being tested, control= the remaining 10 clones). While this could have been automated, this was done manually by copying tables typed in excel into a text file using `vim.` Started building .Rmd of RNA seq analysis. There are issues with the results tables and the MA plot is not made.
