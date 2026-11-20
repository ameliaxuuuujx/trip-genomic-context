#Applying Holm's correction to control for multiple-hypothesis testing. 

#load libraries
library(tidyr)

#Load the bedtools window for each clone. please don't judge the non-looping strategy. 

#####CLONE 16
window_c16 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR16_window.bed", sep = "\t", header = FALSE)
colnames(window_c16) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c16 <- separate(window_c16, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 17
window_c17 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR17_window.bed", sep = "\t", header = FALSE)
colnames(window_c17) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c17 <- separate(window_c17, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 18
window_c18 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR18_window.bed", sep = "\t", header = FALSE)
colnames(window_c18) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c18 <- separate(window_c18, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 19
window_c19 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR19_window.bed", sep = "\t", header = FALSE)
colnames(window_c19) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c19 <- separate(window_c19, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 20
window_c20 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR20_window.bed", sep = "\t", header = FALSE)
colnames(window_c20) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c20 <- separate(window_c20, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 21
window_c21 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR21_window.bed", sep = "\t", header = FALSE)
colnames(window_c21) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c21 <- separate(window_c21, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 22
window_c22 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR22_window.bed", sep = "\t", header = FALSE)
colnames(window_c22) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c22 <- separate(window_c22, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 23
window_c23 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR23_window.bed", sep = "\t", header = FALSE)
colnames(window_c23) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c23 <- separate(window_c23, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 24
window_c24 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR24_window.bed", sep = "\t", header = FALSE)
colnames(window_c24) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c24 <- separate(window_c24, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 25
window_c25 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR25_window.bed", sep = "\t", header = FALSE)
colnames(window_c25) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c25 <- separate(window_c25, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#####CLONE 26
window_c26 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR26_window.bed", sep = "\t", header = FALSE)
colnames(window_c26) <- c("IR_chrom", "IR_start", "IR_end", "barcode", "IR_score", "IR_strand",
                          "gene_chrom", "gene_start", "gene_end", "gene_string", "gene_score", "gene_strand")

window_c26 <- separate(window_c26, gene_string,
                       into = c("gene", "baseMean", "log2FoldChange", "lfcSE",
                                "pvalue", "padj", "symbol", "p_holm"),
                       sep = "_", convert = TRUE)

#extract list of genes within the window of the IRs:
window_list  <- mget(paste0("window_c", 16:26))
window_genes <- unique(unlist(lapply(window_list, `[[`, "gene")))

#length(window_genes) #1288 genes only

#Make a table of p_values from the res_tables in an object called p_values_list

#Load DESeq results table as dataframe
C16 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c16_res_table__df.tsv", sep = "\t")
C17 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c17_res_table__df.tsv", sep = "\t")
C18 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c18_res_table__df.tsv", sep = "\t")
C19 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c19_res_table__df.tsv", sep = "\t")
C20 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c20_res_table__df.tsv", sep = "\t")
C21 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c21_res_table__df.tsv", sep = "\t")
C22 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c22_res_table__df.tsv", sep = "\t")
C23 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c23_res_table__df.tsv", sep = "\t")
C24 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c24_res_table__df.tsv", sep = "\t")
C25 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c25_res_table__df.tsv", sep = "\t")
C26 <- read.delim("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/c26_res_table__df.tsv", sep = "\t")


#Create a list of each clone 
res_list <- list(C16, C17, C18, C19, C20, C21, C22, C23, C24, C25, C26)
#Name the objects in the list 
names(res_list) <- c("C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26")

#Confirm that the gene order is identical for each clone first
all(sapply(res_list, function(res) identical(res[[1]], C16[[1]])))
#[1] TRUE

# extract gene column from one object to label columns by ENSID
all_p_values <- data.frame(gene = C16$gene)

# add pvalue column from each object with the object name as the column header
for (name in names(res_list)) {
  all_p_values[[name]] <- res_list[[name]]$pvalue
}

#Save pvals to compare
pre_holms <- all_p_values

#backup of pre-filtered
pre_filtered_all_p_vals <- all_p_values

##NEW STEP
all_p_values <- all_p_values[all_p_values$gene %in% window_genes, ]

# pool all p-values together
all_p_values_vector <- unlist(all_p_values[, -1])


# apply holm's correction to the whole pool
adjusted_p <- p.adjust(all_p_values_vector, method = "holm")

#Check distrubution of the pvals
#Most of these genes are not differentially expressed after adj testing. I think this is reasonable for RNA seq of TRIP clones

hist(adjusted_p, breaks = 50,
     main = "Distribution of Holm adjusted p-values",
     xlab = "p-value")
sum(adjusted_p < 0.05, na.rm = TRUE)

#For every clone, add a column to the res_table.tsv of the C{clone}_holm. 
adjusted_p_matrix <- matrix(adjusted_p, nrow = nrow(all_p_values), ncol = 11, byrow = FALSE)
adjusted_df <- data.frame(gene = all_p_values$gene, adjusted_p_matrix)
colnames(adjusted_df)[-1] <- paste0(names(res_list), "_holm")

for (name in names(res_list)) {
  res_list[[name]]$p_holm <- adjusted_df[[paste0(name, "_holm")]][
    match(res_list[[name]]$gene, adjusted_df$gene)]
}

# extract back to individual objects
C16 <- res_list[["C16"]]
C17 <- res_list[["C17"]]
C18 <- res_list[["C18"]]
C19 <- res_list[["C19"]]
C20 <- res_list[["C20"]]
C21 <- res_list[["C21"]]
C22 <- res_list[["C22"]]
C23 <- res_list[["C23"]]
C24 <- res_list[["C24"]]
C25 <- res_list[["C25"]]
C26 <- res_list[["C26"]]

#Now that each table has padj and holms_padj I want to compare significance for each clone. 
# check how many significant genes are lost after holm correction

sig_comparison <- data.frame(
  clone = names(res_list),
  sig_padj = sapply(names(res_list), function(name) {
    sum(res_list[[name]]$padj < 0.05, na.rm = TRUE)
  }),
  sig_holm = sapply(names(res_list), function(name) {
    sum(res_list[[name]]$p_holm < 0.05, na.rm = TRUE)
  })
)

# add column for how many were lost
sig_comparison$lost <- sig_comparison$sig_padj - sig_comparison$sig_holm

print(sig_comparison)
#I think this is more promising than try1


#Save the final data tables

#While I can do this not in a loop, this does it instead in one command
for (name in names(res_list)) {
  write.table(res_list[[name]], 
              file = paste0("/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/", 
                            name, "_holms_res_table__df.tsv"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}

#Update bedtools window files
for (cl in 16:26) {
  obj <- paste0("window_c", cl)
  df  <- get(obj)
  df$p_holm <- adjusted_df[[paste0("C", cl, "_holm")]][match(df$gene, adjusted_df$gene)]
  assign(obj, df)
  
  write.table(df,
              file = paste0("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/IR", cl, "_window_holm.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
}



