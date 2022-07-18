library(Seurat)
library(tximport)
library(fishpond)
library(plyr)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(EnsDb.Hsapiens.v79)
#BiocManager::install("EnsDb.Hsapiens.v79")

#load matrix
#without decoys
files<-file.path("/projectnb/bf528/users/hedgehog_2022/project_4/Data_curator/output_salmon/alevin/quants_mat.gz")
txi <- tximport(files, type="alevin")
txicts <- txi$counts 

#convert geneid to symbol
panc <- CreateSeuratObject(counts = txicts, project = "proj4_panc", min.cells = 3, min.features = 200)
# panc@assays$RNA@counts@Dimnames[[1]] <- geneIDs$SYMBOL
# panc@assays$RNA@data@Dimnames[[1]] <- geneIDs$SYMBOL

#---------3.1------------
panc[["percent.mt"]] <- PercentageFeatureSet(panc, pattern = "^MT-")

VlnPlot(panc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(panc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(panc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#filter
panc <- subset(panc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#---------3.2----------
panc <- NormalizeData(panc, normalization.method = "LogNormalize", scale.factor = 10000)
panc <- FindVariableFeatures(panc, selection.method = "vst", nfeatures = 2000)

#top 10 highly variable genes
top10 <- head(VariableFeatures(panc), 10)
gene_names <- sub("[.][0-9]*$", "", top10)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= gene_names, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
top10 <- geneIDs$SYMBOL

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(panc)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot3
plot4
#scaling the data
all.genes <- rownames(panc)
panc <- ScaleData(panc, features = all.genes)

#lineat dim reduction
panc <- RunPCA(panc, features = VariableFeatures(object = panc))
print(panc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(panc, dims = 1:2, reduction = "pca")
DimPlot(panc, reduction = "pca")
DimHeatmap(panc, dims = 1, cells = 500, balanced = TRUE)

#Determine Dimensionality 
panc <- JackStraw(panc, num.replicate = 100)
panc <- ScoreJackStraw(panc, dims = 1:20)
JackStrawPlot(panc, dims=1:20)

#Plot elbow plot
ElbowPlot(panc)

#------------------3.3---------------
#determine clusters
panc <- FindNeighbors(panc, dims = 1:9)
panc <- FindClusters(panc, resolution = 0.5)

#run non-linear dimensional reduction
panc <- RunUMAP(panc, dims=1:10)
DimPlot(panc, reduction = "umap")

#save rds
saveRDS(panc, "panc.rds")

###################################
#Create pie chart to display relative proportions of count in each cluster
#stores the cluster number and count of cell number
cluster_counts <- as.data.frame(table(Idents(panc)))
cluster_counts <- paste(cluster_counts$Var1, cluster_counts$Freq, sep=" - ")

#stores the cluster number and the relative proportions
cluster_proportions <- as.data.frame(prop.table(table(Idents(panc))))
cluster_proportions$Freq <- cluster_proportions$Freq %>% `*`(100) %>% round(2)
cluster_proportions$Percents <- cluster_proportions$Freq
cluster_proportions$Percents <- paste(cluster_proportions$Percents, "%", sep="")

#plot the piechart
cols = colorRampPalette(brewer.pal(8, "Dark2"))(14)
pie(cluster_proportions$Freq, labels=cluster_proportions$Percents, col=cols, main = "Relative Proportions of Cell Numbers\nFor the Identified Clusters")
legend(1.45, 1, legend=cluster_counts, cex = 0.8, fill = cols, bty="n", title="Clusters")



