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
  group = group_vector
)

#filter out minimally expressed genes
num_filtered <- filtered_subset[rowSums(filtered_subset) > 10, ]

#create deseq2 object
dds <- DESeqDataSetFromMatrix(countData = num_filtered,
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
res4 <- as.data.frame(res)[,c("stat","pvalue","padj")]
res_shrink <- lfcShrink(dds, coef="group_control_vs_bipolar", type="apeglm")

#create .csv file
write.csv(as.data.frame(res_shrink), file = "DESeq2_results_shrunk.csv")
read.csv("DESeq2_results_shrunk.csv")

#create averages for conditions
average_counts <- data.frame(
  Bipolar = rowMeans(normalized_counts[, grepl("B", colnames(normalized_counts)), drop = FALSE]),  # Average of Bipolar samples
  Control = rowMeans(normalized_counts[, grepl("C", colnames(normalized_counts)), drop = FALSE])   # Average of Control samples
)

#create heat map
average_counts_s <- scale(average_counts)

View(average_counts_s)

png("heatmap.png", width = 800, height = 600)
pheatmap(average_counts_s,
         cluster_rows = TRUE,          # cluster rows (genes)
         show_rownames = TRUE,         # show gene names
         show_colnames = TRUE,         # show sample names
         fontsize_row = 6,             # adjust row font size
         fontsize_col = 8,             # adjust column font size
         main = "Bipolar vs Control DE Heat Map")
View(heatmap_data_scaled)
dev.off()
dev.flush()
dev.flush()