#Applying Holm's correction to control for multiple-hypothesis testing. 


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
colnames(adjusted_p_matrix) <- paste0(names(res_list), "_holm")

# now add back to res_list
for (name in names(res_list)) {
  res_list[[name]]$p_holm <- adjusted_p_matrix[, paste0(name, "_holm")]
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

#This QC step is alarming because I did not expect more statisticlaly significant genes after multiple hypothesis testing. 

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




