# Changed only the paths


#PJSC20_KM = 4 donors combined, human transcripts only
#In this analysis, I am working on 4 donors combined dataset but we will separate only human transcripts prior normalization.
#The reason is because we saw only few human genes were up-regulated in infected clusters. We reasoned that this might be due to too much
#microsporidian transcripts in the cells causing low expression level of the human genes.



#Load libraries
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

#Processing strategies:
#We will separate human transcripts from each donor.Then normalize the reads prior to integration of 4 donors

#Donor 1:
raw.data <- Read10X(data.dir="/Users/ellamolyukbenasuly/Downloads/GSE268707_RAW/Donor1")
raw.gene.count <- raw.data$`Gene Expression`  #38,612 features, 10,100 cells
raw.HTO.count <- raw.data$`Antibody Capture`  #4 features, 10,100 cells
joint.bcs <- intersect(colnames(raw.gene.count), colnames(raw.HTO.count))
raw.gene.count <- raw.gene.count[,joint.bcs] 
raw.HTO.count <- as.matrix(raw.HTO.count[,joint.bcs])
row.names(raw.HTO.count) <- c('ctrl','1dpi','2dpi','3dpi')
gene.hashtag <- CreateSeuratObject(counts = raw.gene.count)
percent.mt <- PercentageFeatureSet(gene.hashtag, pattern = "^MT-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.hashtag, pattern = "^Eint-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.microsporidia, col.name = "percent.microsporidia")
#separate human transcripts
options(max.print = 100000)
rownames(gene.hashtag)
#Human transcripts = 1:36601
#Microsporidian transcripts = 36602:38612
Donor1.human <- subset(gene.hashtag[1:36601,])
rownames(Donor1.human)
Donor1.human <- NormalizeData(Donor1.human, normalization.method = "LogNormalize", scale.factor = 10000)
Donor1.human <- FindVariableFeatures(Donor1.human, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Donor1.human)
Donor1.human <- ScaleData(Donor1.human, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
Donor1.human[["HTO"]] <- CreateAssayObject(counts=raw.HTO.count)
Donor1.human <- MULTIseqDemux(Donor1.human, assay = "HTO")
table(Donor1.human$MULTI_ID)
#1dpi     2dpi     3dpi     ctrl  Doublet Negative 
#1559     1631     2136     2870      878     1026 
Donor1.human <- subset(Donor1.human, idents = c("ctrl","1dpi","2dpi","3dpi"))
dim(Donor1.human) #8,196 cells
Donor1.human <- subset(Donor1.human, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 2000 & nCount_RNA < 40000 & percent.mt <20)
dim(Donor1.human)
#8,087 cells -> There are 10 cells lesser than processing human+microspordian transcripts together

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)

#Donor 2:
raw.data <- Read10X(data.dir = "/Users/ellamolyukbenasuly/Downloads/GSE268707_RAW/Donor2")
raw.gene.count <- raw.data$`Gene Expression`  #38,612 features and 12,279 cells
raw.HTO.count <- raw.data$`Antibody Capture`  #4 features and 12,279 cells
joint.bcs <- intersect(colnames(raw.gene.count), colnames(raw.HTO.count))
raw.gene.count <- raw.gene.count[,joint.bcs]
raw.HTO.count <- as.matrix(raw.HTO.count[,joint.bcs])
row.names(raw.HTO.count) <- c("ctrl","3hpi","12hpi","1dpi")
gene.hashtag <- CreateSeuratObject(counts = raw.gene.count)
percent.mt <- PercentageFeatureSet(gene.hashtag, pattern = "^MT-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.hashtag, pattern = "^Eint-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.microsporidia, col.name = "percent.microsporidia")
#separate human transcripts
options(max.print = 100000)
rownames(gene.hashtag)
#Human transcripts = 1:36601
#Microsporidian transcripts = 36602:38612
Donor2.human <- subset(gene.hashtag[1:36601,])
rownames(Donor2.human)
Donor2.human <- NormalizeData(Donor2.human, normalization.method = "LogNormalize", scale.factor = 10000)
Donor2.human <- FindVariableFeatures(Donor2.human, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Donor2.human)
Donor2.human <- ScaleData(Donor2.human, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
Donor2.human[["HTO"]] <- CreateAssayObject(counts = raw.HTO.count)
Donor2.human <- MULTIseqDemux(Donor2.human, assay = "HTO")
table(Donor2.human$MULTI_ID)
#12hpi     1dpi     3hpi     ctrl  Doublet Negative 
#3483     3219     3434      415      628     1100
remove(gene.hashtag)
Donor2.human <- subset(Donor2.human, idents = c("ctrl","3hpi","12hpi","1dpi"))
Donor2.human <- subset(Donor2.human, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 2200 & nCount_RNA < 35000 & percent.mt <20)
dim(Donor2.human)
#10,013 cells -> 49 cells less than combined transcripts

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)

#Donor 3 and Donor 4
raw.data <- Read10X(data.dir="/Users/ellamolyukbenasuly/Downloads/GSE268707_RAW/Donor34")
raw.gene.count <- raw.data$`Gene Expression`  #38,612 features and 30,355 cells
raw.HTO.count <- raw.data$`Antibody Capture`  #6 features and 30,355 cells
joint.bcs <- intersect(colnames(raw.gene.count), colnames(raw.HTO.count))
raw.gene.count <- raw.gene.count[,joint.bcs]
raw.HTO.count <- as.matrix(raw.HTO.count[,joint.bcs])
rownames(raw.HTO.count) <- c("D3_ctrl","D3_1dpi","D3_2dpi","D4_ctrl","D4_1dpi","D4_2dpi")
gene.hashtag <- CreateSeuratObject(counts = raw.gene.count)
percent.mt <- PercentageFeatureSet(gene.hashtag, pattern = "^MT-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.hashtag, pattern = "^Eint-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.microsporidia, col.name = "percent.microsporidia")
#separate human transcripts
options(max.print = 100000)
rownames(gene.hashtag)
#Human transcripts = 1:36601
#Microsporidian transcripts = 36602:38612
Donor3_4.human <- subset(gene.hashtag[1:36601,])
rownames(Donor3_4.human)
Donor3_4.human <- NormalizeData(Donor3_4.human, normalization.method = "LogNormalize", scale.factor = 10000)
Donor3_4.human <- FindVariableFeatures(Donor3_4.human, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Donor3_4.human)
Donor3_4.human <- ScaleData(Donor3_4.human, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
Donor3_4.human[["HTO"]] <- CreateAssayObject(counts = raw.HTO.count)
Donor3_4.human <- MULTIseqDemux(Donor3_4.human, assay = "HTO")
table(Donor3_4.human$MULTI_ID)
#D3-1dpi  D3-2dpi  D3-ctrl  D4-1dpi  D4-2dpi  D4-ctrl  Doublet Negative 
#3500     5984      830     5583     4431     1572     6535     1920 

Donor3.human <- subset(Donor3_4.human, idents = c("D3-ctrl","D3-1dpi","D3-2dpi"))
Donor4.human <- subset(Donor3_4.human, idents = c("D4-ctrl","D4-1dpi","D4-2dpi"))

Donor3.human <- subset(Donor3.human, subset =  nFeature_RNA > 400 & nFeature_RNA < 4000 & nCount_RNA > 1000 & nCount_RNA < 20000 & percent.mt <20)
dim(Donor3.human) 
#10,206 cells remain from 10,254 cells of combined dataset

Donor4.human <- subset(Donor4.human, subset =  nFeature_RNA > 400 & nFeature_RNA <3800 & nCount_RNA > 1000 & nCount_RNA < 15000 & percent.mt <20)
dim(Donor4.human)
#11,427 cells remain from 11,468 cells of combined dataset

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)
remove(Donor3_4.human)

#Assign donor name into each cell and store the data in 'Exp' slot
Donor1.human$Exp <- "Donor1"
Donor2.human$Exp <- "Donor2"
Donor3.human$Exp <- "Donor3"
Donor4.human$Exp <- "Donor4"

#Perform integration and batch correction
#First, select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(Donor1.human, Donor2.human, Donor3.human, Donor4.human), nfeatures = 3000)
options(max.print = 100000)
features

#Identify anchors which will be used for integration
combined.anchors <- FindIntegrationAnchors(object.list = list(Donor1.human, Donor2.human, Donor3.human, Donor4.human), anchor.features = features)

remove(Donor1.human)
remove(Donor2.human)
remove(Donor3.human)
remove(Donor4.human)

#Combined datasets
combined.dataset <- IntegrateData(anchorset = combined.anchors)
dim(combined.dataset)
#39,733 cells with 3,000 features from 39,881 of combined dataset

#Next, let's see the integration results
DefaultAssay(combined.dataset) <- "integrated"
combined.dataset <- ScaleData(combined.dataset, verbose = FALSE)
combined.dataset <- RunPCA(combined.dataset)
percent.pc <- combined.dataset@reductions$pca@stdev / sum(combined.dataset@reductions$pca@stdev) * 100
cum.percent.pc <- cumsum(percent.pc)
cum.percent.pc
#[1]   6.270388  11.294934  16.220203  20.707770  24.725748  28.112483  30.727803  33.211383
#[9]  35.609745  37.914314  40.141590  42.304463  44.357861  46.355328  48.325793  50.210062
#[17]  52.069857  53.833117  55.579620  57.260624  58.909239  60.547362  62.162328  63.739551
#[25]  65.298021  66.844840  68.380343  69.899667  71.405720  72.896115  74.352516  75.795610
#[33]  77.212149  78.625565  80.030782  81.429446  82.818673  84.195373  85.559847  86.911473
#[41]  88.259771  89.600768  90.933117  92.248584  93.560822  94.858603  96.156175  97.443019
#[49]  98.722512 100.000000
#PC1-PC43

#Run UMAP
combined.dataset <- RunUMAP(combined.dataset, dims = 1:43)
Idents(combined.dataset) <- "Exp"

png("PJSC20_UMAP_all_donors.png", width = 12, height = 10, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", pt.size = 0.1, cols = "Set1")
dev.off()

png("PJSC20_UMAP_all_donors_split.png", width = 40, height = 12, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", pt.size = 0.1, split.by = "Exp", ncol = 4, cols = "Set1")
dev.off()

#Let's see where are the infected cells
png("PJSC15_UMAP_percent.microsporidia.png", width = 12, height = 10, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = "percent.microsporidia")
dev.off()
#Look like the bottom cluster contains the most infected cells. This is similar to C4 that we found in D1 and D2

#Perform the clustering
set.seed(2022)
Idents(combined.dataset) <- "orig.ident"
combined.dataset <- FindNeighbors(combined.dataset, dims = 1:43)

#previously tested 0.5 to 0.1 resolutions and chose to move forward with 0.2 

#Resolution 0.2
combined.dataset <- FindClusters(combined.dataset, resolution = 0.2)
png("PJSC20_Dimplot_0.2_FigureQ_NoLabel.png", width = 12, height = 10, units = "cm", res = 300)
DimPlot(combined.dataset, pt.size = 0.1)
dev.off()

png("PJSC20_Dimplot_0.2_Control_2.png", width = 50, height = 15, units = "cm", res = 300)
DimPlot(combined.dataset, pt.size = 0.1, split.by = "MULTI_ID")
dev.off()


combined.dataset$MULTI_ID <- recode_factor(combined.dataset$MULTI_ID, "D3-ctrl" = "ctrl", "D3-1dpi" = "1dpi", "D3-2dpi" = "2dpi", "D4-ctrl"  = "ctrl", "D4-1dpi" = "1dpi", "D4-2dpi" = "2dpi")
combined.dataset$MULTI_ID <- factor(combined.dataset$MULTI_ID, levels = c("ctrl","3hpi","12hpi","1dpi","2dpi","3dpi"))

png("PJSC20_UMAP_all_donors_split_0.2.png", width = 17, height = 20, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "MULTI_ID", pt.size = 0.1) + facet_wrap(~combined.dataset$MULTI_ID + combined.dataset$Exp)
dev.off()
#9 clusters (1 singletons identified)

png("PJSC20_UMAP_0.2_percent.microsporidia.png", width = 12, height = 10, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = "percent.microsporidia", label = T)
dev.off()

png("PJSC20_UMAP_split_Timepoint_NoLabel.png", width = 20, height = 10, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "MULTI_ID", pt.size = 0.1)
dev.off ()

#Quantification of infected cells
combined.dataset <- FindClusters(combined.dataset, resolution = 0.20)
cells_with_microsporidia_2 <- subset(combined.dataset, subset = percent.microsporidia > 2)
infected_cell_barcodes_2 <- colnames(cells_with_microsporidia_2)
infected_cell_barcodes_2 <- data.frame(infected_cell_barcodes_2)
dim(infected_cell_barcodes_2) #6,908 cells (148 cells reduced from the combined dataset)
combined.dataset$infected_cell_2 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_2$infected_cell_barcodes_2, "Yes", "No" )

write.csv(infected_cell_barcodes_2, "PJSC20_infected_cell_barcodes_2.csv")


infected_cell_barcodes_and_clusters <- combined.dataset@meta.data[,c(5,10,11,13)]
PJSC20_infected_cell_barcodes_and_clusters_DF <- data.frame(infected_cell_barcodes_and_clusters)

write.csv(PJSC20_infected_cell_barcodes_and_clusters_DF, "PJSC20_infected_cell_barcodes_and_clusters_DF.csv")

Idents(combined.dataset) <- "infected_cell_2"
png("Infected_cell_0.2.png", width = 12, height = 10, units = "cm", res = 300)
DimPlot(combined.dataset, pt.size = 0.1, label.size = 5, cols = c("gray","mediumblue"))
dev.off()

png("Quanification_Infected_cell_0.2.png", width = 12, height = 10, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.2, infected_cell_2) %>% count() %>% group_by(integrated_snn_res.0.2) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.2, y=percent, fill = infected_cell_2)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 17, face = "plain")) + scale_fill_brewer(palette = "Greys")
dev.off()

#Export infected cells barcodes and their human clusters
view(combined.dataset)
PJSC20_infected_cell_barcodes_and_clusters <- combined.dataset@meta.data[,c(5,10,11,13)]
write.csv(combined.dataset@meta.data, "PJSC20_infected_cell_barcodes_and_clusters.csv")

Macrophage.signature <- c("FABP4","PPARG","APOE","APOC1","DDX60","SPP1","CSF3R","IL1B","RGS2","ISG15","CHIT1","STAB1","TOP2A","WFDC2")
png("PJS20_DotPlot_MO_signatures_blue.png", width = 13, height = 10, units = "cm", res = 300)
DotPlot(combined.dataset, features = Macrophage.signature,cols = c("white","blue")) + RotatedAxis() + coord_flip()
dev.off()

install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scRNAseq")

library(SingleR)
library(celldex)
library(scRNAseq)

Human.primary.cells <- HumanPrimaryCellAtlasData()
Human.primary.cells

dataset.for.singleR <- GetAssayData(combined.dataset)
dim(dataset.for.singleR)

Predicted.dataset.main <- SingleR(test = dataset.for.singleR, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.main)
Predicted.dataset.main

table(Predicted.dataset.main$labels)

Predicted.dataset <- SingleR(test = dataset.for.singleR, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
Predicted.dataset
table(Predicted.dataset$labels)

plotScoreHeatmap(Predicted.dataset)

DefaultAssay(combined.dataset) <- "RNA"
png("Macrophage_sigs_Feature_Plot_2.png", width = 30, height = 25, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = c("FABP4","PPARG","APOE","APOC1","DDX60","SPP1","CSF3R","IL1B","RGS2","ISG15","CHIT1","STAB1","TOP2A","WFDC2"), ncol = 4, label = T, label.size = 5)
dev.off()

subset_C0 <- subset(combined.dataset, subset = integrated_snn_res.0.2 == 0)
dim(subset_C0)
singleR_C0 <- GetAssayData(subset_C0)
Predicted.dataset.C0 <- SingleR(test = singleR_C0, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
table(Predicted.dataset.C0$labels)
write.csv(table(Predicted.dataset.C0$labels), "Predicted.dataset.C0.csv")

subset_C1 <- subset(combined.dataset, subset = integrated_snn_res.0.2 == 1)
dim(subset_C1)
singleR_C1 <- GetAssayData(subset_C1)
Predicted.dataset.C1 <- SingleR(test = singleR_C1, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
table(Predicted.dataset.C1$labels)
write.csv(table(Predicted.dataset.C1$labels), "Predicted.dataset.C1.csv")

subset_C2 <- subset(combined.dataset, subset = integrated_snn_res.0.2 == 2)
dim(subset_C2)
singleR_C2 <- GetAssayData(subset_C2)
Predicted.dataset.C2 <- SingleR(test = singleR_C2, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
table(Predicted.dataset.C2$labels)
write.csv(table(Predicted.dataset.C2$labels), "Predicted.dataset.C2.csv")

subset_C3 <- subset(combined.dataset, subset = integrated_snn_res.0.2 == 3)
dim(subset_C3)
singleR_C3 <- GetAssayData(subset_C3)
Predicted.dataset.C3 <- SingleR(test = singleR_C3, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
table(Predicted.dataset.C3$labels)
write.csv(table(Predicted.dataset.C3$labels), "Predicted.dataset.C3.csv")

subset_C4 <- subset(combined.dataset, subset = integrated_snn_res.0.2 == 4)
dim(subset_C4)
singleR_C4 <- GetAssayData(subset_C4)
Predicted.dataset.C4 <- SingleR(test = singleR_C4, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
table(Predicted.dataset.C4$labels)
write.csv(table(Predicted.dataset.C4$labels), "Predicted.dataset.C4.csv")

subset_C5 <- subset(combined.dataset, subset = integrated_snn_res.0.2 == 5)
dim(subset_C5)
singleR_C5 <- GetAssayData(subset_C5)
Predicted.dataset.C5 <- SingleR(test = singleR_C5, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
table(Predicted.dataset.C5$labels)
write.csv(table(Predicted.dataset.C5$labels), "Predicted.dataset.C5.csv")

subset_C6 <- subset(combined.dataset, subset = integrated_snn_res.0.2 == 6)
dim(subset_C6)
singleR_C6 <- GetAssayData(subset_C6)
Predicted.dataset.C6 <- SingleR(test = singleR_C6, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
table(Predicted.dataset.C6$labels)
write.csv(table(Predicted.dataset.C6$labels), "Predicted.dataset.C6.csv")

subset_C7 <- subset(combined.dataset, subset = integrated_snn_res.0.2 == 7)
dim(subset_C7)
singleR_C7 <- GetAssayData(subset_C7)
Predicted.dataset.C7 <- SingleR(test = singleR_C7, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
table(Predicted.dataset.C7$labels)
write.csv(table(Predicted.dataset.C7$labels), "Predicted.dataset.C7.csv")

subset_C8 <- subset(combined.dataset, subset = integrated_snn_res.0.2 == 8)
dim(subset_C8)
singleR_C8 <- GetAssayData(subset_C8)
Predicted.dataset.C8 <- SingleR(test = singleR_C8, ref = Human.primary.cells, assay.type.test = 1, labels = Human.primary.cells$label.fine)
table(Predicted.dataset.C8$labels)
write.csv(table(Predicted.dataset.C8$labels), "Predicted.dataset.C8.csv")

png("Human_Cell_Atlas_grey.png", width = 18, height = 12, units = "cm", res = 300)
ggplot(All_Clusters_Final, aes(X, Y, fill= Percent)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="gray8") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1.5))
dev.off()

png("Human_Cell_Atlas_100323.png", width = 18, height = 12, units = "cm", res = 300)
ggplot(All_Clusters_Final_10032023, aes(X, Y, fill= Percent)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1.5))
dev.off()

DefaultAssay(combined.dataset) <- "RNA"
RNA_allmarkers <- FindAllMarkers(combined.dataset, min.pct = 0.25)
write.csv(RNA_allmarkers, "RNA_allmarkers.csv")

write.csv(combined.dataset@meta.data, "PJSC20_MetaData.csv")

RNA_allmarkers.top10 <- RNA_allmarkers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
png("PJSC20_DotPlot_RNA_top10_Flip_0.2_rdbu.png", width = 50, height = 11, units = "cm", res = 300)
DotPlot(combined.dataset, features = unique(RNA_allmarkers.top10$gene), cols = "RdBu") + RotatedAxis()
dev.off()


C2_vs_C0 <- FindMarkers(combined.dataset, ident.1 = "2", ident.2 = "0")
write.csv(C2_vs_C0, "PJSC20_C2_vs_C0.csv")

C2_vs_C1 <- FindMarkers(combined.dataset, ident.1 = "2", ident.2 = "1")
write.csv(C2_vs_C1, "PJSC20_C2_vs_C1.csv")

C0_vs_C1 <- FindMarkers(combined.dataset, ident.1 = "0", ident.2 = "1")
write.csv(C0_vs_C1, "PJSC20_C0_vs_C1.csv")

C3_vs_C0 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "0")
write.csv(C3_vs_C0, "PJSC20_C3_vs_C0.csv")

C3_vs_C1 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "1")
write.csv(C3_vs_C1, "PJSC20_C3_vs_C1.csv")

C3_vs_C1 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "1")
write.csv(C3_vs_C1, "PJSC20_C3_vs_C1.csv")

C3_vs_C2 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "2")
write.csv(C3_vs_C2, "PJSC20_C3_vs_C2.csv")

C3_vs_C4 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "4")
write.csv(C3_vs_C4, "PJSC20_C3_vs_C4.csv")

C3_vs_C5 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "5")
write.csv(C3_vs_C5, "PJSC20_C3_vs_C5.csv")

C3_vs_C6 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "6")
write.csv(C3_vs_C6, "PJSC20_C3_vs_C6.csv")

C3_vs_C7 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "7")
write.csv(C3_vs_C7, "PJSC20_C3_vs_C7.csv")

C3_vs_C8 <- FindMarkers(combined.dataset, ident.1 = "3", ident.2 = "8")
write.csv(C3_vs_C8, "PJSC20_C3_vs_C8.csv")

C4_vs_C0 <- FindMarkers(combined.dataset, ident.1 = "4", ident.2 = "0")
write.csv(C4_vs_C0, "PJSC20_C4_vs_C0.csv")

C4_vs_C1 <- FindMarkers(combined.dataset, ident.1 = "4", ident.2 = "1")
write.csv(C4_vs_C1, "PJSC20_C4_vs_C1.csv")

C4_vs_C2 <- FindMarkers(combined.dataset, ident.1 = "4", ident.2 = "2")
write.csv(C4_vs_C2, "PJSC20_C4_vs_C2.csv")

C4_vs_C3 <- FindMarkers(combined.dataset, ident.1 = "4", ident.2 = "3")
write.csv(C4_vs_C3, "PJSC20_C4_vs_C3.csv")

C4_vs_C5 <- FindMarkers(combined.dataset, ident.1 = "4", ident.2 = "5")
write.csv(C4_vs_C5, "PJSC20_C4_vs_C5.csv")

C4_vs_C6 <- FindMarkers(combined.dataset, ident.1 = "4", ident.2 = "6")
write.csv(C4_vs_C6, "PJSC20_C4_vs_C6.csv")

C4_vs_C7 <- FindMarkers(combined.dataset, ident.1 = "4", ident.2 = "7")
write.csv(C4_vs_C7, "PJSC20_C4_vs_C7.csv")

C4_vs_C8 <- FindMarkers(combined.dataset, ident.1 = "4", ident.2 = "8")
write.csv(C4_vs_C8, "PJSC20_C4_vs_C8.csv")

C5_vs_C0 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "0")
write.csv(C5_vs_C0, "PJSC20_C5_vs_C0.csv")

C5_vs_C1 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "1")
write.csv(C5_vs_C1, "PJSC20_C5_vs_C1.csv")

C5_vs_C2 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "2")
write.csv(C5_vs_C2, "PJSC20_C5_vs_C2.csv")

C5_vs_C3 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "3")
write.csv(C5_vs_C3, "PJSC20_C5_vs_C3.csv")

C5_vs_C4 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "4")
write.csv(C5_vs_C4, "PJSC20_C5_vs_C4.csv")

C5_vs_C6 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "6")
write.csv(C5_vs_C6, "PJSC20_C5_vs_C6.csv")

C5_vs_C7 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "7")
write.csv(C5_vs_C7, "PJSC20_C5_vs_C7.csv")

C5_vs_C8 <- FindMarkers(combined.dataset, ident.1 = "5", ident.2 = "8")
write.csv(C5_vs_C8, "PJSC20_C5_vs_C8.csv")

png("Macrophage_sigs_Feature_Plot.png", width = 40, height = 20, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = c("TNF", "IL10"), ncol = 5, label = T, label.size = 5)
dev.off()

png("PJSC20_quantification_timepoint_in_each_cluster.png", width = 16, height = 10, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.2, MULTI_ID) %>% count() %>% group_by(integrated_snn_res.0.2) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.2, y=percent, fill = MULTI_ID)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 15, face = "plain")) + scale_fill_brewer(palette = "Set1")
dev.off()

png("PJSC20_Res0.20_quantification_donor.png", width = 14, height = 8, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.2, Exp) %>% count() %>% group_by(integrated_snn_res.0.2) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.2, y=percent, fill = Exp)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 15, face = "plain")) + scale_fill_brewer(palette = "Set1")
dev.off()

Macrophage.signature_2 <- c("FABP4","INHBA","PPARG","APOE","APOC1","DDX60","SPP1","MERTK","LGMN","SIGLEC10","CD14","INSIG1","OSM","THBS1","MRC1","FCN1","MARCO","CSF3R","IL1B","RGS2","ISG15","CHIT1","STAB1","TOP2A","WFDC2","NFKB1","STAT1","STAT2","TFEB","NR1H3","PPARA","CREB1","CEBPB","CCL2","CCL3","CXCL10","AM2","GPR183","CCL13","TREM2","TGFB1")
png("PJS20_DotPlot_MO_signatures_Updated.png", width = 13, height = 10, units = "cm", res = 300)
DotPlot(combined.dataset, features = Macrophage.signature_2,cols = c("white","blue")) + RotatedAxis() + coord_flip()
dev.off()

Cyto_Chemokines <- c("CXCL10", "CXCL11","CCL5","CCR7","IDO7","TGFB","CCL14","CCL22","SRB1","PPARG")
png("PJS20_DotPlot_MO_signatures_Updated.png", width = 13, height = 10, units = "cm", res = 300)
DotPlot(combined.dataset, features =Cyto_Chemokines,cols = c("white","blue")) + RotatedAxis() + coord_flip()
dev.off()

Immune_reg <- c("GPR183","CCL13","TREM2","TGFB1","SPP1")
DotPlot(combined.dataset, features = Immune_reg,cols = c("white","blue")) + RotatedAxis() + coord_flip()
dev.off()

Mac_Sigs <- c("FCN1","FABP4","INHBA","SSP1","MERTK","LGMN","SIGLEC10")
DotPlot(combined.dataset, features = Immune_reg,cols = c("white","blue")) + RotatedAxis() + coord_flip()
dev.off()

Signatures_nature_paper <- c("CDKN2A","RND3","SULT1B1","CTSF","APOE","APOC1","EPHX1","SDC2","SPP1","CENPF","DDX60","FCGR38","MARCO","S100A8","S100A9","CCL23","IL1R2","CSF3R","IL1B","STMN1","CCL2","CXCR4","IL4I1","CHIT1","HLA-DRB1","VM01","CHD9","ADA2","KIF1C","HMGB2")
DotPlot(combined.dataset, features = Signatures_nature_paper,cols = c("white","blue")) + RotatedAxis() + coord_flip()
dev.off()

combined.dataset@meta.data %>% group_by(MULTI_ID, infected_cell_2) %>% count() %>% group_by(MULTI_ID) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=MULTI_ID, y=percent, fill = infected_cell_2)) + geom_col() + 
  ylab("Percent") + xlab("TimePoints") + theme(text = element_text(size = 17, face = "plain")) + scale_fill_brewer(palette = "Set1") + ylim(0, 50)

DefaultAssay(combined.dataset) <- "RNA"
png("PJSC20_VlnPlot_C3andC5_markers.png", width = 25, height = 9, units = "cm", res = 300)
VlnPlot(combined.dataset, features = c("ISG15","MX1","IFIT1","IFIT3","RSAD2","INHBA","IER3","CIR1","MALAT1","NEAT1"),ncol = 5, pt.size = 0) & theme(axis.title.x = element_blank())
dev.off()
######################################
infected_timepoints_donor <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_2, combined.dataset@meta.data$Exp)
infected_timepoints_donor
write.csv(infected_timepoints_donor,"infected_timepoints_donor.csv")
#Donor 1
#       No  Yes
#ctrl  2825   17
#3hpi     0    0
#12hpi    0    0
#1dpi  1224  310
#2dpi  1087  511
#3dpi  1656  457

#Donor 2
#       No  Yes
#ctrl   395    2
#3hpi  3238   45
#12hpi 2802  468
#1dpi  2470  593
#2dpi     0    0
#3dpi     0    0

#Donor 3
#       No  Yes
#ctrl   815   11
#3hpi     0    0
#12hpi    0    0
#1dpi  2690  771
#2dpi  4776 1143
#3dpi     0    0

#Donor 4
#       No  Yes
#ctrl  1539   13
#3hpi     0    0
#12hpi    0    0
#1dpi  3989 1530
#2dpi  3319 1037
#3dpi     0    0

no.cells.timepoints <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$orig.ident, combined.dataset@meta.data$Exp)
no.cells.timepoints
#Donor 1
#ctrl           2842
#3hpi              0
#12hpi             0
#1dpi           1534
#2dpi           1598
#3dpi           2113

#Donor 2
#ctrl            397
#3hpi           3283
#12hpi          3270
#1dpi           3063
#2dpi              0
#3dpi              0

#Donor 3
#ctrl            826
#3hpi              0
#12hpi             0
#1dpi           3461
#2dpi           5919
#3dpi              0

#Donor 4
#ctrl           1552
#3hpi              0
#12hpi             0
#1dpi           5519
#2dpi           4356
#3dpi              0

percent_infected <- combined.dataset@meta.data %>% group_by(integrated_snn_res.0.2, infected_cell_2) %>% count() %>% group_by(integrated_snn_res.0.2) %>% mutate(percent=100*n/sum(n)) 

write.csv(percent_infected, "percent_infected_table.csv")

no.cells.clusters <- table(combined.dataset@meta.data$integrated_snn_res.0.2, combined.dataset@meta.data$orig.ident)
no.cells.clusters
#    SeuratProject
#0          10817
#1          10496
#2          10214
#3           4948
#4           1791
#5           1068
#6            182
#7            168
#8             76

microsporidian_transcripts_table <- combined.dataset@meta.data %>% group_by(integrated_snn_res.0.2, percent.microsporidia) %>% count() %>% group_by(integrated_snn_res.0.2) %>% mutate(percent=100*n/sum(n)) 
write.csv(microsporidian_transcripts_table, "microsporidian_transcripts_table.csv")

mitochondria_table <- combined.dataset@meta.data %>% group_by(integrated_snn_res.0.2, percent.mt) %>% count() %>% group_by(integrated_snn_res.0.2) %>% mutate(percent=100*n/sum(n)) 
write.csv(mitochondria_table, "mitochondria_table.csv")

nFeature_table <- combined.dataset@meta.data %>% group_by(integrated_snn_res.0.2, nFeature_RNA) %>% count() %>% group_by(integrated_snn_res.0.2) %>% mutate(percent=100*n/sum(n)) 
write.csv(nFeature_table, "nFeature_table.csv")

p1 <- ggplot(combined.dataset@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + geom_point() + scale_color_gradient(low = "grey90", high = "black") + stat_smooth(method = Im) + theme_classic() + ylab("nGenes") + xlab("nUMIs") + theme(text = element_text(size = 15, face = "bold"))
p1

p2 <- ggplot(combined.dataset@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.microsporidia)) + geom_point() + scale_color_gradient(low = "grey90", high = "black") + stat_smooth(method = Im) + theme_classic() + ylab("nGenes") + xlab("nUMIs") + theme(text = element_text(size = 15, face = "bold"))
p2

write.csv(combined.dataset@meta.data, "combined_dataset_metadata_table.csv")

install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
BiocManager::install("SeuratWrappers")
library(monocle3)
library(cicero)
library(SeuratWrappers)

alldata.cds <- as.cell_data_set(combined.dataset)

DefaultAssay(combined.dataset) <- "RNA"
png("Interferons.png", width = 70, height = 50, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = c("IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNA10","IFNA13","IFNA14","IFNA17","IFNA21","IFNB1","IFNG","IFNK","IFNL1","IFNL2","IFNL3","IFNW1","IL6","IL29","IL28A","IL28B","IRF3","IRF7","IRF1","IRF5","IRF8","IRF9","TNF","NFKB1","JAK1","JAK2","TYK2","STAT1","STAT2"), ncol = 7, label = T, label.size = 5)
dev.off()

infected_timepoints <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_2, combined.dataset@meta.data$Exp, combined.dataset@meta.data$integrated_snn_res.0.2)
infected_timepoints
write.csv(infected_timepoints,"infected_timepoints_per_cluster_donor.csv")

no.total.cells.clusters <- table(combined.dataset@meta.data$integrated_snn_res.0.2, combined.dataset@meta.data$orig.ident,combined.dataset@meta.data$Exp,combined.dataset@meta.data$MULTI_ID)
no.total.cells.clusters
write.csv(no.total.cells.clusters,"no.total.cells.clusters.csv")


FeaturePlot(combined.dataset, features = c("TLR1","TLR2","TLR3","TLR4","TLR5","TLR6","TLR7","TLR8","TLR9","TLR10"), label = T, label.size = 5)

FeaturePlot(combined.dataset, features = c("TNFR1", "TNFR2"), label = T, label.size = 5)
