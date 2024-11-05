#install DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#install tidyr & dplyr
install.packages("tidyverse")

#install pheatmap
install.packages("pheatmap")

#load libraries
library(DESeq2)
library(tidyr)
library(dplyr)
library(pheatmap)

#load gene table
filtered <- read.table("C:/Users/ashlu/OneDrive/Desktop/bioinfo final project/GSE42546_CleanedRawCounts.txt.gz",
  header = TRUE, sep = "\t", row.names = 1)
columns_to_keep <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10",
                     "B11", "B12", "B13", "B14", "B15", "B16", 
                     "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
                     "C11", "C12", "C13", "C14", "C15", "C16")
  filtered_subset <- filtered %>%
    select(all_of(columns_to_keep))

#create table of samples x group
  x <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10",
         "B11", "B12", "B13", "B14", "B15", "B16", 
         "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
         "C11", "C12", "C13", "C14", "C15", "C16")
  y <- c(rep("bipolar", 16), rep("control", 16))
df <- data.frame(sample = x, group = y)

#define group vector
group_vector <- factor(y)

#create dataset object
dataset_object <- data.frame(
  sample = colnames(filtered_subset),
  group = group_vector)

#create deseq2 object
dds <- DESeqDataSetFromMatrix(countData = filtered_subset,
    colData = dataset_object,
    design=~group)
  print(dds)
  
#normalize data
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)

#results
coef_names <- resultsNames(dds)
print(coef_names)
res <- results(dds, name="group_control_vs_bipolar")
res_shrink <- lfcShrink(dds, coef="group_control_vs_bipolar", type="apeglm")

#create .csv file
write.csv(as.data.frame(res_shrink), file = "DESeq2_results_shrunk.csv")
read.csv("DESeq2_results_shrunk.csv")

#create heat map
gene_ids <- rownames(normalized_counts)
heatmap_data <- normalized_counts[gene_ids, , drop = FALSE]
heatmap_data_scaled <- t(scale(t(heatmap_data)))
if (any(is.na(heatmap_data_scaled)) || any(is.infinite(heatmap_data_scaled))) {
  stop("Scaled data contains NA or infinite values. Please check the input data.")
}
png("heatmap.png", width = 800, height = 600)
pheatmap(heatmap_data_scaled,
         cluster_rows = TRUE,          # Cluster rows (genes)
         cluster_cols = TRUE,          # Cluster columns (samples)
         show_rownames = TRUE,         # Show gene names
         show_colnames = TRUE,         # Show sample names
         fontsize_row = 6,             # Adjust row font size
         fontsize_col = 8,             # Adjust column font size
         main = "Bipolar vs Control - All Samples Heat Map")
dev.off()