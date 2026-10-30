#Make bedfiles for the RNA seq files and for the IR clones
library(rtracklayer)
library(dplyr)
library(tidyverse)
#GTF: https://www.gencodegenes.org/mouse/release_M25.html 

#Make bed files for each clone IR locations
C16_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C16_holms_res_table__df.tsv', sep = "\t")
C17_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C17_holms_res_table__df.tsv', sep = "\t")
C18_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C18_holms_res_table__df.tsv', sep = "\t")
C19_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C19_holms_res_table__df.tsv', sep = "\t")
C20_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C20_holms_res_table__df.tsv', sep = "\t")
C21_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C21_holms_res_table__df.tsv', sep = "\t")
C22_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C22_holms_res_table__df.tsv', sep = "\t")
C23_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C23_holms_res_table__df.tsv', sep = "\t")
C24_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C24_holms_res_table__df.tsv', sep = "\t")
C25_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C25_holms_res_table__df.tsv', sep = "\t")
C26_genes <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/results/compile_deseq/C26_holms_res_table__df.tsv', sep = "\t")

#Load GTF file 
GTF <- import('/Users/cmdb/trip-genomic-context/rnaseq/data/gencode.vM25.annotation.gtf')
GTF_df <- as.data.frame(GTF)

#extract gene ID and bed file
GTF_genes <- GTF_df %>%
  filter(type == "gene") %>%
  dplyr::select(
    seqnames,    # chromosome number
    start,
    end,
    strand,
    gene_id, # ENSID
    gene_type, 
    gene_name
  ) %>%
  # Clean ENSID to match your table
  # Remove version number from ENSID
  mutate(gene_id_clean = gsub("\\..*", "", gene_id))

########LOOP 
# Define your clone list
clone_list <- paste0("C", seq(16, 26))

# Empty summary table
summary_table <- data.frame(
  clone = character(),
  total_genes = numeric(),
  genes_with_coordinates = numeric(),
  genes_lost = numeric(),
  percent_lost = numeric(),
  stringsAsFactors = FALSE
)

# Loop through clones
for (clone in clone_list) {
  
  #Loop over the results tables named {clone}_genes
  genes_dt <- get(paste0(clone, "_genes"))
  
  #Add gene coordinates 
  merged <- genes_dt %>%
    left_join(GTF_genes,
              by = c("gene" = "gene_id_clean"))
  
  # Create BED
  bed <- merged %>%
    mutate(
      name = paste(gene,
                   baseMean,
                   log2FoldChange,
                   lfcSE,
                   pvalue,
                   padj,
                   symbol,
                   p_holm,
                   sep = "_"),
      score = 0,
      start = start - 1
    ) %>%
    dplyr::select(
      seqnames,
      start,
      end,
      name,
      score,
      strand
    )
  
  # Count lost genes
  total <- nrow(genes_dt)
  with_coords <- sum(!is.na(bed$seqnames))
  lost <- sum(is.na(bed$seqnames))
  
  # Save bed object dynamically
  assign(paste0(clone, "_bed"), bed)
  
  # Save bed file
  bed_clean <- bed %>%
    filter(!is.na(seqnames))
  
  write.table(bed_clean,
              file = paste0("/Users/cmdb/trip-genomic-context/rnaseq/results/bedfiles/",
                            clone, "genes.bed"),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
  
  # Add to summary table
  summary_table <- rbind(summary_table,
                         data.frame(
                           clone = clone,
                           total_genes = total,
                           genes_with_coordinates = with_coords,
                           genes_lost = lost,
                           percent_lost = round(lost / total * 100, 1)
                         ))
  
  message(paste(clone, "done -", lost, "genes lost"))
}

# View summary
print(summary_table)

#################################

#IR clones
trip_data <- read.delim('/Users/cmdb/trip-genomic-context/rnaseq/data/results/bedfiles/trip_clone_IRs.txt',
                        header = TRUE,
                        sep = "\t")

# Extract clone number from TRIP_cell_line column
trip_data <- trip_data %>%
  mutate(clone = gsub("EStTA-PB-B-", "", TRIP_cell_line))

# list the clones
clones <- unique(trip_data$clone)

for (clone in clones) {
  clone_data <- trip_data %>%
    filter(clone == !!clone) %>%
    mutate(
      start = position - 1,
      end = position,
      name = barcode,
      score = 0
    ) %>%
    dplyr::select(
      chromosome,    # chr
      start,         # start 0-based
      end,           # end
      name,          # barcode as name
      score,         # normalization counts as score
      strand         # strand
    )
  # Save as IR{clone}.bed
  write.table(clone_data,
              file = paste0("/Users/cmdb/trip-genomic-context/rnaseq/results/IR_", 
                            clone, ".bed"),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE)
}
