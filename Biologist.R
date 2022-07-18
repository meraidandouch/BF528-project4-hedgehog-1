# Very simple simple script made for creating a txt file with a list
# of genes in each cluster that meet the p-val threshold.

library(tidyverse)

# Read in list of genes
genes <- read_csv("cell_markers.csv")

# Set threshold
threshold <- 1e-20
# Filter by adjusted p-value
filtered_genes <- filter(genes, p_val_adj <= threshold)
# Filtering by adjusted p-value is problematic as the number of cells
# in each cluster that meet this parameter vary pretty drstically

# Create vector of clusters
clusters <- unique(filtered_genes$cluster)
# For each cluster output a text file with the names of the genes
for (cluster1 in clusters) {
  temp <- filter(filtered_genes, cluster == cluster1)
  write.table(temp$gene, paste('cluster', cluster1, '.txt', sep=""),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}