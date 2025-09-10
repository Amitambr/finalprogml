# Changed only the paths
#Load necessary libraries

library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

#Donor 1:
raw.data <- Read10X(data.dir="/Users/ellamolyukbenasuly/Downloads/GSE268707_RAW/Donor1")
raw.gene.count <- raw.data$`Gene Expression`  #38,612 features, 10,100 cells
raw.HTO.count <- raw.data$`Antibody Capture`  #4 features, 10,100 cells
joint.bcs <- intersect(colnames(raw.gene.count), colnames(raw.HTO.count))
raw.gene.count <- raw.gene.count[,joint.bcs] 
raw.HTO.count <- as.matrix(raw.HTO.count[,joint.bcs])
row.names(raw.HTO.count) <- c('ctrl','1dpi','2dpi','3dpi')
gene.hashtag <- CreateSeuratObject(counts = raw.gene.count)
gene.hashtag <- NormalizeData(gene.hashtag, normalization.method = "LogNormalize", scale.factor = 10000)
gene.hashtag <- FindVariableFeatures(gene.hashtag, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(gene.hashtag)
gene.hashtag <- ScaleData(gene.hashtag, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
gene.hashtag[["HTO"]] <- CreateAssayObject(counts=raw.HTO.count)
gene.hashtag <- MULTIseqDemux(gene.hashtag, assay = "HTO")
table(gene.hashtag$MULTI_ID)
#1dpi     2dpi     3dpi     ctrl  Doublet Negative 
#1559     1631     2136     2870      878     1026 
gene.singlet <- subset(gene.hashtag, idents = c("ctrl","1dpi","2dpi","3dpi"))
dim(gene.singlet) #8,196 cells
percent.mt <- PercentageFeatureSet(gene.singlet, pattern = "^MT-")
gene.singlet <- AddMetaData(gene.singlet, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.singlet, pattern = "^Eint-")
gene.singlet <- AddMetaData(gene.singlet, metadata = percent.microsporidia, col.name = "percent.microsporidia")
Donor1 <- subset(gene.singlet, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 2000 & nCount_RNA < 40000 & percent.mt <20)
dim(Donor1)
#8,097 cells remain

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
gene.hashtag <- NormalizeData(gene.hashtag, normalization.method = "LogNormalize", scale.factor = 10000)
gene.hashtag <- FindVariableFeatures(gene.hashtag, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(gene.hashtag)
gene.hashtag <- ScaleData(gene.hashtag, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
gene.hashtag[["HTO"]] <- CreateAssayObject(counts = raw.HTO.count)
gene.hashtag <- MULTIseqDemux(gene.hashtag, assay = "HTO")
table(gene.hashtag$MULTI_ID)
#12hpi     1dpi     3hpi     ctrl  Doublet Negative 
#3483     3219     3434      415      628     1100
gene.singlet <- subset(gene.hashtag, idents = c("ctrl","3hpi","12hpi","1dpi"))
percent.mt <- PercentageFeatureSet(gene.singlet, pattern = "^MT-")
gene.singlet <- AddMetaData(gene.singlet, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.singlet, pattern = "^Eint-")
gene.singlet <- AddMetaData(gene.singlet, metadata = percent.microsporidia, col.name = "percent.microsporidia")
Donor2 <- subset(gene.singlet, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & nCount_RNA > 2200 & nCount_RNA < 35000 & percent.mt <20)
dim(Donor2)
#10,062 cells remain

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)

#Save Donor 1 and Donor 2 as RDS
saveRDS(Donor1, "Donor1.rds")
saveRDS(Donor2, "Donor2.rds")

#Donor 3 and Donor 4:
raw.data <- Read10X(data.dir="/Users/ellamolyukbenasuly/Downloads/GSE268707_RAW/Donor34")
raw.gene.count <- raw.data$`Gene Expression`  #38,612 features and 30,355 cells
raw.HTO.count <- raw.data$`Antibody Capture`  #6 features and 30,355 cells
joint.bcs <- intersect(colnames(raw.gene.count), colnames(raw.HTO.count))
raw.gene.count <- raw.gene.count[,joint.bcs]
raw.HTO.count <- as.matrix(raw.HTO.count[,joint.bcs])
rownames(raw.HTO.count) <- c("D3_ctrl","D3_1dpi","D3_2dpi","D4_ctrl","D4_1dpi","D4_2dpi")
gene.hashtag <- CreateSeuratObject(counts = raw.gene.count)
gene.hashtag <- NormalizeData(gene.hashtag, normalization.method = "LogNormalize", scale.factor = 10000)
gene.hashtag <- FindVariableFeatures(gene.hashtag, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(gene.hashtag)
gene.hashtag <- ScaleData(gene.hashtag, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
gene.hashtag[["HTO"]] <- CreateAssayObject(counts = raw.HTO.count)
gene.hashtag <- MULTIseqDemux(gene.hashtag, assay = "HTO")
table(gene.hashtag$MULTI_ID)
#D3-1dpi  D3-2dpi  D3-ctrl  D4-1dpi  D4-2dpi  D4-ctrl  Doublet Negative 
#3500     5984      830     5583     4431     1572     6535     1920 
percent.mt <- PercentageFeatureSet(gene.hashtag, pattern = "^MT-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.hashtag, pattern = "^Eint-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.microsporidia, col.name = "percent.microsporidia")

gene.singlet3 <- subset(gene.hashtag, idents = c("D3-ctrl","D3-1dpi","D3-2dpi"))
gene.singlet4 <- subset(gene.hashtag, idents = c("D4-ctrl","D4-1dpi","D4-2dpi"))

Donor3 <- subset(gene.singlet3, subset =  nFeature_RNA > 400 & nFeature_RNA < 4000 & nCount_RNA > 1000 & nCount_RNA < 20000 & percent.mt <20)
dim(Donor3) 
#10,254 cells remain

Donor4 <- subset(gene.singlet4, subset =  nFeature_RNA > 400 & nFeature_RNA <3800 & nCount_RNA > 1000 & nCount_RNA < 15000 & percent.mt <20)
dim(Donor4)
#11,468 cells remain

#Save Donor 3 and Donor 4 as RDS
saveRDS(Donor3, "Donor3.rds")
saveRDS(Donor4, "Donor4.rds")

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet3)
remove(gene.singlet4)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)

#Assign donor name into each cell and store the data in 'Exp' slot
Donor1$Exp <- "Donor1"
Donor2$Exp <- "Donor2"
Donor3$Exp <- "Donor3"
Donor4$Exp <- "Donor4"


#Next, we are going to perform integration and batch correction
#First, select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(Donor1, Donor2, Donor3, Donor4), nfeatures = 3000)
options(max.print = 100000)
features

#Identify the anchors
#We then identify anchors using the FindIntegrationAnchors() function, which will take a list of Seurat objects as input, and use these anchors to integrate the four datasets together with IntegrateData().
combined.anchors <- FindIntegrationAnchors(object.list = list(Donor1, Donor2, Donor3, Donor4), anchor.features = features)
#CCA is performed to find an anchor
#23,047 anchors are found -> 4,206 anchors remain after filtering
#It will take ~1 hr to run

#Warning message:
#In CheckDuplicateCellNames(object.list = object.list) :
#Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.

remove(Donor1)
remove(Donor2)
remove(Donor3)
remove(Donor4)

#Merging 4 datasets
combined.dataset <- IntegrateData(anchorset = combined.anchors)
dim(combined.dataset)
#3,000 features, 39,881 cells

#Next, let's see the integration results
DefaultAssay(combined.dataset) <- "integrated"
combined.dataset <- ScaleData(combined.dataset, verbose = FALSE)
combined.dataset <- RunPCA(combined.dataset)
percent.pc <- combined.dataset@reductions$pca@stdev / sum(combined.dataset@reductions$pca@stdev) * 100
cum.percent.pc <- cumsum(percent.pc)
cum.percent.pc
#[1]  21.96502  26.90431  31.45381  34.76272  38.01465  40.94041  43.64577  46.26728  48.33232  50.32710  52.16523
#[12]  53.87988  55.55547  57.17248  58.70376  60.21769  61.71437  63.17951  64.62870  66.03014  67.37471  68.69862
#[23]  69.97698  71.22500  72.45574  73.67421  74.88090  76.05985  77.22107  78.37429  79.51928  80.65682  81.78651
#[34]  82.90236  84.01715  85.12363  86.21831  87.30997  88.39688  89.47771  90.54750  91.61526  92.67512  93.73334
#[45]  94.78968  95.84098  96.88502  97.92733  98.96453 100.00000
#PC1-PC41


#Run UMAP
combined.dataset <- RunUMAP(combined.dataset, dims = 1:41)
Idents(combined.dataset) <- "Exp"

png("PJSC22_UMAP_all_donors.png", width = 12, height = 10, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", pt.size = 0.1, cols = "Set1")
dev.off()

png("PJSC22_UMAP_all_donors_split.png", width = 24, height = 20, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", pt.size = 0.1, split.by = "Exp", ncol = 2, cols = "Set1")
dev.off()

#Split according to the infection timepoints
#First, we need to change the name of MULTI_ID for D3 and D4 from "D3-ctrl" to onlt "ctrl"
combined.dataset$MULTI_ID <- recode_factor(combined.dataset$MULTI_ID, "D3-ctrl"  = "ctrl", "D3-1dpi" = "1dpi", "D3-2dpi" = "2dpi", "D4-ctrl"  = "ctrl", "D4-1dpi" = "1dpi", "D4-2dpi" = "2dpi")

#Then, order the tinepoints
combined.dataset$MULTI_ID <- factor(combined.dataset$MULTI_ID, levels = c("ctrl","3hpi","12hpi","1dpi","2dpi","3dpi"))

png("PJSC14_UMAP_all_donors_split2.png", width = 28, height = 20, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "MULTI_ID", cols = "Set1", pt.size = 0.1) + facet_wrap(~combined.dataset$MULTI_ID + combined.dataset$Exp)
dev.off()

png("PJSC14_UMAP_Timpoint_Label.png", width = 20, height = 15, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", split.by = "MULTI_ID", ncol = 3, label = T, pt.size = 0.1) 
dev.off()

#Find number of cells in each timepoints
no.cells.timepoints <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$orig.ident, combined.dataset@meta.data$Exp)
no.cells.timepoints


#Perform the clustering
set.seed(2022)
Idents(combined.dataset) <- "orig.ident"
combined.dataset <- FindNeighbors(combined.dataset, dims = 1:41)

#Resolution 0.4
combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
png("PJSC22_UMAP_res0.4.png", width = 12, height = 10, units = "cm", res = 300)
DimPlot(combined.dataset, pt.size = 0.1, label = T, label.size = 4)
dev.off()
#12 clusters

png("PJSC22_UMAP_res0.4_byDonor_Nolabel.png", width = 20, height = 20, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", split.by = "Exp", ncol = 2, label = F, pt.size = 0.1)
dev.off()

DimPlot(combined.dataset, reduction = "umap", split.by = "MULTI_ID", ncol = 2, label = F, pt.size = 0.1)

png("PJSC22_UMAP_all_donors_split_0.4.png", width = 17, height = 20, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "MULTI_ID", pt.size = 0.1) + facet_wrap(~combined.dataset$MULTI_ID + combined.dataset$Exp, ncol = 4, label = T)
dev.off()

png("PJSC20_UMAP_0.4_percent.microsporidia.png", width = 12, height = 10, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = "percent.microsporidia", label = T)
dev.off()

#Quantification of infected cells
combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_2 <- subset(combined.dataset, subset = percent.microsporidia > 2)
infected_cell_barcodes_2 <- colnames(cells_with_microsporidia_2)
infected_cell_barcodes_2 <- data.frame(infected_cell_barcodes_2)
dim(infected_cell_barcodes_2) #7056 cells
combined.dataset$infected_cell_2 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_2$infected_cell_barcodes_2, "Yes", "No" )


PJSC22_infected_cell_barcodes_and_clusters <- combined.dataset@meta.data[,c(9,10,13,14)]
write.csv(PJSC22_infected_cell_barcodes_and_clusters, "PJSC22_infected_cell_barcodes_and_clusters.csv")

pwcIdents(combined.dataset) <- "infected_cell_2"
png("Infected_cell_0.4_byDonor.png", width = 20, height = 20, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "Exp", ncol = 2, pt.size = 0.1, cols = c("gray","mediumblue"))
dev.off()

png("Quanification_Infected_cell_0.4.png", width = 15, height = 10, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.4, infected_cell_2) %>% count() %>% group_by(integrated_snn_res.0.4) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.4, y=percent, fill = infected_cell_2)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 17, face = "plain")) + scale_fill_brewer(palette = "Greys")
dev.off()

png("PJSC22_Res0.4_quantification.png", width = 15, height = 8, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.4, Exp) %>% count() %>% group_by(integrated_snn_res.0.4) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.4, y=percent, fill = Exp)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 15, face = "plain")) + scale_fill_brewer(palette = "Set1")
dev.off()

#I am going forward with Resolution 0.4

#Let's see in each cluster where are the cells from
png("PJSC14_Res0.4_quantification.png", width = 14, height = 8, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.4, Exp) %>% count() %>% group_by(integrated_snn_res.0.4) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.4, y=percent, fill = Exp)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 15, face = "plain")) + scale_fill_manual(values=c("darkorange2","deepskyblue4","darkolivegreen4","brown"))
dev.off()
#cluster 7,9,11 are C3 and C4 specific

#Which cluster contains microsporidian transcripts
Idents(combined.dataset) <- "integrated_snn_res.0.4"
png("PJSC14_Res0.4_microsporidian.png", width = 12, height = 10, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = "percent.microsporidia", label = T, pt.size = 0.1)
dev.off()

#Quantification of infected cells 
cells_with_microsporidia_2 <- subset(combined.dataset, subset = percent.microsporidia > 2)
infected_cell_barcodes_2 <- colnames(cells_with_microsporidia_2)
infected_cell_barcodes_2 <- data.frame(infected_cell_barcodes_2)
dim(infected_cell_barcodes_2) #7,056 cells
combined.dataset$infected_cell_2 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_2$infected_cell_barcodes_2, "Yes", "No" )

Idents(combined.dataset) <- "integrated_snn_res.0.4"
p1 <- DimPlot(combined.dataset, label = T, cols = "Paired", pt.size = 0.1, label.size = 5)

Idents(combined.dataset) <- "infected_cell_2"
p2 <- DimPlot(combined.dataset, pt.size = 0.1, label.size = 5, cols = c("gray","mediumblue"))

p3 <- combined.dataset@meta.data %>% group_by(integrated_snn_res.0.4, infected_cell_2) %>% count() %>% group_by(integrated_snn_res.0.4) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.4, y=percent, fill = infected_cell_2)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 17, face = "plain")) + scale_fill_brewer(palette = "Greys")

png("PJSC14_Res0.4_microsporidian_quantification.png", width = 38, height = 9, units = "cm", res = 300)
p1+p2+p3
dev.off()

png("PJSC22_Res0.4_quantification_Timepoint.png", width = 15, height = 8, units = "cm", res = 300)
combined.dataset@meta.data %>% group_by(integrated_snn_res.0.4, MULTI_ID) %>% count() %>% group_by(integrated_snn_res.0.4) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.4, y=percent, fill = MULTI_ID)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 15, face = "plain")) + scale_fill_brewer(palette = "Set1")
dev.off()

#Look at top10 markers of each cluster
allmarkers.top10 <- allmarkers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

#Plot heatmap of these top10 markers
png("PJSC14_Heatmap_top10.png", width = 38, height = 20, units = "cm", res = 300)
DoHeatmap(combined.dataset, features = unique(allmarkers.top10$gene))
dev.off()

#top5markers
allmarkers.top5 <- allmarkers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
png("PJSC14_Heatmap_top5.png", width = 10, height = 5, units = "cm", res = 300)
DoHeatmap(combined.dataset, features = unique(allmarkers.top5$gene))
dev.off()

FeaturePlot(combined.dataset, features = unique(allmarkers.top5$gene), label = T)

#From D1 and D2 data processing, we found D2 specific clusters. These clusters are highly expressed in IFN induced genes such as
#CCL8, CCL13, ISG15, IFIT2, CXCL10
#I would like to see where these genes are expressed when we combined 4 donors
DefaultAssay(combined.dataset) <- "RNA"
png("PJSC22_Feature_C3.png", width = 18, height = 20, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = c("CCL8", "IFIT3","ISG15", "IFIT2", "CXCL10"), label = T, pt.size = 0.1)
dev.off()

png("PJSC22_Violin_C3.png", width = 18, height = 20, units = "cm", res = 300)
VlnPlot(combined.dataset, features = c("CCL8", "IFIT3","ISG15", "IFIT2", "CXCL10"))
dev.off()

#I looked at the tutorial from Seurat https://satijalab.org/seurat/articles/integration_introduction.html
#It suggests switching to "RNA" slot prior doing DGE analysis. Previously, I did the DGE on "integrated" slot.
DefaultAssay(combined.dataset) <- "RNA"
RNA_allmarkers <- FindAllMarkers(combined.dataset, min.pct = 0.25)
RNA_allmarkers.top5 <- RNA_allmarkers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

write.csv(RNA_allmarkers,"PJSC22_RNA_allmarker.csv")

write.csv(combined.dataset@meta.data,"PJSC22_MetaData.csv")

png("PJSC22_DotPlot_RNA_top5_vertical_RdBu_35_25.png", width = 25, height = 35, units = "cm", res = 300)
DotPlot(combined.dataset, features = unique(RNA_allmarkers.top5$gene), cols = "RdBu") + RotatedAxis() + coord_flip()
dev.off()

RNA_allmarkers.top10 <- RNA_allmarkers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
png("PJSC14_DotPlot_RNA_top10.png", width = 20, height = 40, units = "cm", res = 300)
DotPlot(combined.dataset, features = unique(RNA_allmarkers.top10$gene)) + RotatedAxis()
dev.off()

DefaultAssay(combined.dataset) <- "integrated"
DoHeatmap(combined.dataset, features = unique(RNA_allmarkers.top5$gene)) + scale_fill_gradientn(colours = rev(brewer.pal(n=9, name="RdBu")), na.value = "white")

#2022-02-16
#From these integration, we could identify 2 mega-clusters, namely one without microsporidian transcript (Right), and one with microsporidian transcripts (Left)
#To answer how host cells response to microsporidia infection we look at the host genes in infected clusters (C4, C5, C6, C9)
#When I look at the markers for C4, most of the host genes in C4 are down-regulated 
#No host gene that are up-regulated. This might be because we normalized both microsporidian and human transcripts together. 
#Hence, we need to separate between human and microsporidian transcripts for the data-processing



#We found that Cluster 7 and Cluster 9 are Donor 3 and Donor 4 specific clusters and they have similar expression profiles
#But Cluster 9 contains microsporidian transcripts but cluster 7 does not
#Now let's see what are the gene markers for C7 and C9


#dDoHeatmap(combined.dataset, features = unique(RNA_allmarkers.top5$gene))

#How many microsporidian genes are expressed in C4,5,6 and 9
Idents(combined.dataset) <- "integrated_snn_res.0.4"

combined.dataset@meta.data %>% group_by(integrated_snn_res.0.4, combined.dataset$) %>% count() %>% group_by(integrated_snn_res.0.4) %>% mutate(percent=100*n/sum(n)) %>% ungroup() %>% ggplot(aes(x=integrated_snn_res.0.4, y=percent, fill = infected_cell_2)) + geom_col() + ylab("Percent") + xlab("Clusters") + theme(text = element_text(size = 17, face = "plain")) + scale_fill_brewer(palette = "Greys")

#######################
#2022-05-22
#Plot dot plots only infected clusters (C4,5,6,9)
infected_cluster_markers <- RNA_allmarkers.top5$cluster == c("4","5","6","9")
infected_cluster_markers <- RNA_allmarkers.top5[c(41:50,51:60,61:70,91:100),]

png("PJSC14_DotPlot_RNA_top10_infected_clusters.png", width = 43, height = 6, units = "cm", res = 300)
DotPlot(combined.dataset, features = unique(infected_cluster_markers$gene), idents = c("4","5","6","9")) + RotatedAxis()
dev.off()

DefaultAssay(combined.dataset) <- "integrated"
DoHeatmap(combined.dataset, features = infected_cluster_markers$gene)

DoHeatmap(combined.dataset, features = infected_cluster_markers$gene)

infected_timepoints <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_2, combined.dataset@meta.data$Exp)
infected_timepoints
write.csv(infected_timepoints,"infected_timepoints.csv")
#Donor 1
#       No  Yes
#ctrl  2825   17
#3hpi     0    0
#12hpi    0    0
#1dpi  1224  312
#2dpi  1087  517
#3dpi  1656  459

#Donor 2
#       No  Yes
#ctrl   395    2
#3hpi  3238   45
#12hpi 2802  479
#1dpi  2470  631
#2dpi     0    0
#3dpi     0    0

#Donor 3
#       No  Yes
#ctrl   815   10
#3hpi     0    0
#12hpi    0    0
#1dpi  2690  775
#2dpi  4776 1188
#3dpi     0    0

#Donor 4
#       No  Yes
#ctrl  1539   10
#3hpi     0    0
#12hpi    0    0
#1dpi  3989 1524
#2dpi  3319 1087
#3dpi     0    0

no.cells.timepoints <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$orig.ident, combined.dataset@meta.data$Exp)
no.cells.timepoints
#Donor 1
#ctrl           2842
#3hpi              0
#12hpi             0
#1dpi           1536
#2dpi           1604
#3dpi           2115

#Donor 2
#ctrl            397
#3hpi           3283
#12hpi          3281
#1dpi           3101
#2dpi              0
#3dpi              0

#Donor 3
#ctrl            825
#3hpi              0
#12hpi             0
#1dpi           3465
#2dpi           5964
#3dpi              0

#Donor 4
#ctrl           1549
#3hpi              0
#12hpi             0
#1dpi           5513
#2dpi           4406
#3dpi              0

percent_infected <- combined.dataset@meta.data %>% group_by(integrated_snn_res.0.4, infected_cell_2) %>% count() %>% group_by(integrated_snn_res.0.4) %>% mutate(percent=100*n/sum(n)) 

write.csv(percent_infected, "percent_infected_table.csv")

no.cells.clusters <- table(combined.dataset@meta.data$integrated_snn_res.0.4, combined.dataset@meta.data$orig.ident)
no.cells.clusters
#    SeuratProject
#0          10844
#1           8096
#2           6433
#3           4197
#4           2908
#5           2814
#6           1547
#7           1368
#8           1298
#9            181
#10           133
#11            62

no.cells.timepoints <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$orig.ident, combined.dataset@meta.data$Exp,combined.dataset@meta.data$integrated_snn_res.0.4)
no.cells.timepoints
write.csv(no.cells.timepoints,"Donors_Timepoint_Cluster.csv")

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_1 <- subset(combined.dataset, subset = percent.microsporidia > 1)
infected_cell_barcodes_1 <- colnames(cells_with_microsporidia_1)
infected_cell_barcodes_1 <- data.frame(infected_cell_barcodes_1)
dim(infected_cell_barcodes_1) #7797 cells
combined.dataset$infected_cell_1 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_1$infected_cell_barcodes_1, "Yes", "No" )

infected_timepoints_1 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_1)
infected_timepoints_1

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_2 <- subset(combined.dataset, subset = percent.microsporidia > 2)
infected_cell_barcodes_2 <- colnames(cells_with_microsporidia_2)
infected_cell_barcodes_2 <- data.frame(infected_cell_barcodes_2)
dim(infected_cell_barcodes_2) #7056 cells
combined.dataset$infected_cell_2 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_2$infected_cell_barcodes_2, "Yes", "No" )

infected_timepoints_2 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_2)
infected_timepoints_2

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_3 <- subset(combined.dataset, subset = percent.microsporidia > 3)
infected_cell_barcodes_3 <- colnames(cells_with_microsporidia_3)
infected_cell_barcodes_3 <- data.frame(infected_cell_barcodes_3)
dim(infected_cell_barcodes_3) #7056 cells
combined.dataset$infected_cell_3 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_3$infected_cell_barcodes_3, "Yes", "No" )

infected_timepoints_3 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_3)
infected_timepoints_3

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_4 <- subset(combined.dataset, subset = percent.microsporidia > 4)
infected_cell_barcodes_4 <- colnames(cells_with_microsporidia_4)
infected_cell_barcodes_4 <- data.frame(infected_cell_barcodes_4)
dim(infected_cell_barcodes_4) #7056 cells
combined.dataset$infected_cell_4 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_4$infected_cell_barcodes_4, "Yes", "No" )

infected_timepoints_4 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_4)
infected_timepoints_4

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_5 <- subset(combined.dataset, subset = percent.microsporidia > 5)
infected_cell_barcodes_5 <- colnames(cells_with_microsporidia_5)
infected_cell_barcodes_5 <- data.frame(infected_cell_barcodes_5)
dim(infected_cell_barcodes_5) #7056 cells
combined.dataset$infected_cell_5 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_5$infected_cell_barcodes_5, "Yes", "No" )

infected_timepoints_5 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_5)
infected_timepoints_5

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_10 <- subset(combined.dataset, subset = percent.microsporidia > 10)
infected_cell_barcodes_10 <- colnames(cells_with_microsporidia_10)
infected_cell_barcodes_10 <- data.frame(infected_cell_barcodes_10)
dim(infected_cell_barcodes_10) #7056 cells
combined.dataset$infected_cell_10 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_10$infected_cell_barcodes_10, "Yes", "No" )

infected_timepoints_10 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_10)
infected_timepoints_10

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_20 <- subset(combined.dataset, subset = percent.microsporidia > 20)
infected_cell_barcodes_20 <- colnames(cells_with_microsporidia_20)
infected_cell_barcodes_20 <- data.frame(infected_cell_barcodes_20)
dim(infected_cell_barcodes_20) #7056 cells
combined.dataset$infected_cell_20 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_20$infected_cell_barcodes_20, "Yes", "No" )

infected_timepoints_20 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_20)
infected_timepoints_20

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_40 <- subset(combined.dataset, subset = percent.microsporidia > 40)
infected_cell_barcodes_40 <- colnames(cells_with_microsporidia_40)
infected_cell_barcodes_40 <- data.frame(infected_cell_barcodes_40)
dim(infected_cell_barcodes_40) #7056 cells
combined.dataset$infected_cell_40 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_40$infected_cell_barcodes_40, "Yes", "No" )

infected_timepoints_40 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_40)
infected_timepoints_40

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_80 <- subset(combined.dataset, subset = percent.microsporidia > 80)
infected_cell_barcodes_80 <- colnames(cells_with_microsporidia_80)
infected_cell_barcodes_80 <- data.frame(infected_cell_barcodes_80)
dim(infected_cell_barcodes_80) #7056 cells
combined.dataset$infected_cell_80 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_80$infected_cell_barcodes_80, "Yes", "No" )

infected_timepoints_80 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_80)
infected_timepoints_80

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_0 <- subset(combined.dataset, subset = percent.microsporidia >= 0)
infected_cell_barcodes_0 <- colnames(cells_with_microsporidia_0)
infected_cell_barcodes_0 <- data.frame(infected_cell_barcodes_0)
dim(infected_cell_barcodes_0) #7056 cells
combined.dataset$infected_cell_0 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_0$infected_cell_barcodes_0, "Yes", "No" )

infected_timepoints_0 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_0)
infected_timepoints_0

combined.dataset <- FindClusters(combined.dataset, resolution = 0.4)
cells_with_microsporidia_0.1 <- subset(combined.dataset, subset = percent.microsporidia >= 0.1)
infected_cell_barcodes_0.1 <- colnames(cells_with_microsporidia_0.1)
infected_cell_barcodes_0.1 <- data.frame(infected_cell_barcodes_0.1)
dim(infected_cell_barcodes_0.1) #10185 cells
combined.dataset$infected_cell_0.1 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_0.1$infected_cell_barcodes_0.1, "Yes", "No" )

infected_timepoints_0.1 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_0.1)
infected_timepoints_0.1

cells_with_microsporidia_0.25 <- subset(combined.dataset, subset = percent.microsporidia >= 0.25)
infected_cell_barcodes_0.25 <- colnames(cells_with_microsporidia_0.25)
infected_cell_barcodes_0.25 <- data.frame(infected_cell_barcodes_0.25)
dim(infected_cell_barcodes_0.25) #8641 cells
combined.dataset$infected_cell_0.25 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_0.25$infected_cell_barcodes_0.25, "Yes", "No" )

infected_timepoints_0.25 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_0.25)
infected_timepoints_0.25

cells_with_microsporidia_0.5 <- subset(combined.dataset, subset = percent.microsporidia >= 0.5)
infected_cell_barcodes_0.5 <- colnames(cells_with_microsporidia_0.5)
infected_cell_barcodes_0.5 <- data.frame(infected_cell_barcodes_0.5)
dim(infected_cell_barcodes_0.5) #8310 cells
combined.dataset$infected_cell_0.5 <- ifelse(rownames(combined.dataset@meta.data) %in% infected_cell_barcodes_0.5$infected_cell_barcodes_0.5, "Yes", "No" )

infected_timepoints_0.5 <- table(combined.dataset@meta.data$MULTI_ID, combined.dataset@meta.data$infected_cell_0.5)
infected_timepoints_0.5
