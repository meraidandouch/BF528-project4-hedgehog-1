library(tidyverse)
library(Seurat)

#pbmc <- readRDS(file = "/projectnb/bf528/project_4_scrnaseq/GSM2230760_seurat.rda")

pbmc <- readRDS("GSM2230760_seurat.rda")

cell_markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)

write.csv(cell_markers,"cell_markers.csv")

top10 <- cell_markers %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt = avg_log2FC) %>% 
  filter(p_val < 0.00001)
write.csv(top10,"top10_cell_markers.csv")


cell_types <- c("Delta", "Beta", "Pancreatic stellate cells", "Acinar", 
                "Alpha", "Ductal", "Beta", "Beta", "Alpha",
                "Pancreatic stellate cells", "Ductal", "Acinar", "Mast")

names(cell_types) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, cell_types)

write.csv(cell_markers,"cell_markers.csv")


cluster_info <- top10 %>% # I used this to find cell_types 
  select(cluster, gene) %>% 
  group_by(cluster) %>%
  mutate(all_genes = paste(gene, collapse = ",")) %>% 
  select(!gene) %>% unique()

# Clustered umap
png("tsne_celltype.png", height = 700, width = 900)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

# Clustered Heatmap
png("heatmap_celltype.png", height = 900, width = 1000)
DoHeatmap(pbmc, features = top10$gene)
dev.off()

top1 <- cell_markers %>% 
  group_by(cluster) %>% 
  top_n(n=1, wt = avg_log2FC) %>% 
  filter(p_val_adj < 0.05)

cluster_top1_info <- top1 %>% # I used this to find cell_types 
  select(cluster, gene) %>% 
  group_by(cluster) %>%
  mutate(all_genes = paste(gene, collapse = ",")) %>% 
  select(!gene) %>% unique()

# gamma 
png("cluster0.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("SST","SP100","PRRG3","LDB1","DNMT3A","RAB12","COL6A3","S1PR2","AC010359.3","SNX27"))
dev.off()

# beta
png("cluster1.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("INS","DLK1","CDKN1C","IAPP","DEPP1","PCSK1","G6PC2","ADCYAP1","TSC22D1","HADH"))
dev.off()

# psc
png("cluster2.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("FN1","COL1A2","BGN","COL1A1","SPARC","C11orf96","COL6A1","SERPINE1","TIMP1","COL4A2"))
dev.off()

# acinar
png("cluster3.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("REG1B","REG1A","CTRB2","PRSS2","SPINK1","CPA1","CELA3A","PRSS1","REG3A","CPB1"))
dev.off()

# alpha
png("cluster4.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("TTR","GCG","TM4SF4","CHGB","ALDH1A1","CLU","IRX2","MUC13","GC","PEG10"))
dev.off()

# ductal
png("cluster5.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("KRT19","MMP7","PMEPA1","ANXA2","TACSTD2","CXCL1","SERPING1","CXCL8","CXCL6","H19"))
dev.off()

# delta
png("cluster6.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("EEF1A2","EDN3","RPS3A","PCSK1N","CITED2","PCSK1","RPL21","RPS25","INS","GAPDH"))
dev.off()

# gamma2
png("cluster7.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("ACER3","RPS6KA5","HELLPAR","ZNF320","SEMA3E","ADAM10","AL022322.2","CHST12","SLC35E3","NABP1"))
dev.off()

# alpha2
png("cluster8.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("GC","PLCE1","TM4SF4","CRYBA2","IGFBP2","GPX3","SLC7A2","TTR","CHGB","CFC1"))
dev.off()

# psc2
png("cluster9.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("COL1A1","COL3A1","COL6A2","COL1A2","SPARC","SERPINE1","BGN","COL6A1","IL11","FN1"))
dev.off()

# ductal2
png("cluster10.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("CRP","KRT18","LCN2","KRT7","CFB","C3","ANXA4","CFTR","KRT8","AL049839.2"))
dev.off()

# acinar
png("cluster11.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("ALDOB","CELA3A","CTRB1","PRSS1","C3","G0S2","CTRB2","SPINK1","PRSS2","AL049839.2"))
dev.off()

png("cluster12.png", height = 2000, width = 2000)
VlnPlot(pbmc, features = c("ACP5","APOE","HLA-DRA","SDS","LAPTM5","ALOX5AP","IFI30","AC007192.1","LYZ","CTSB"))
dev.off()
