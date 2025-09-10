#changed the paths and a little fix marked as FIX

#Load necessary libraries
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

#Processing strategies: we will separate E. intestinalis transcripts from each donor.Then normalize the reads prior to integration of 4 donors

#Donor 1:
raw.data <- Read10X(data.dir="/Users/ellamolyukbenasuly/Downloads/GSE268707_RAW/Donor1")
raw.gene.count <- raw.data$`Gene Expression`  #38,612 features, 10,100 cells
raw.HTO.count <- raw.data$`Antibody Capture` #4 features, 10,100 cells
joint.bcs <- intersect(colnames(raw.gene.count), colnames(raw.HTO.count))
raw.gene.count <- raw.gene.count[,joint.bcs] 
raw.HTO.count <- as.matrix(raw.HTO.count[,joint.bcs])
row.names(raw.HTO.count) <- c('ctrl','1dpi','2dpi','3dpi')
gene.hashtag <- CreateSeuratObject(counts = raw.gene.count)
percent.mt <- PercentageFeatureSet(gene.hashtag, pattern = "^MT-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.mt, col.name = "percent.mt")
percent.microsporidia <- PercentageFeatureSet(gene.hashtag, pattern = "^Eint-")
gene.hashtag <- AddMetaData(gene.hashtag, metadata = percent.microsporidia, col.name = "percent.microsporidia")

#separate E. intestinalis transcripts
options(max.print = 100000)
rownames(gene.hashtag)
#Human transcripts = 1:36601
#Microsporidian transcripts = 36602:38612
Donor1.microsporidian <- subset(gene.hashtag[36602:38612,])
rownames(Donor1.microsporidian)
Donor1.microsporidian <- NormalizeData(Donor1.microsporidian, normalization.method = "LogNormalize", scale.factor = 10000)
Donor1.microsporidian <- FindVariableFeatures(Donor1.microsporidian, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Donor1.microsporidian)
Donor1.microsporidian <- ScaleData(Donor1.microsporidian, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
Donor1.microsporidian[["HTO"]] <- CreateAssayObject(counts=raw.HTO.count)
Donor1.microsporidian <- MULTIseqDemux(Donor1.microsporidian, assay = "HTO")
table(Donor1.microsporidian$MULTI_ID)
#1dpi     2dpi     3dpi     ctrl  Doublet Negative 
#1559     1631     2136     2870      878     1026 
Donor1.microsporidian <- subset(Donor1.microsporidian, idents = c("ctrl","1dpi","2dpi","3dpi"))
dim(Donor1.microsporidian) #8,196 cells

VlnPlot(Donor1.microsporidian, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.microsporidia"))

#### >>> FIX: subset infected cells directly (no external barcode table needed)
Donor1.microsporidian <- subset(Donor1.microsporidian, subset = percent.microsporidia > 2)
dim(Donor1.microsporidian)  # should drop to ~1–2k cells
#### <<<

#### (old code, now unnecessary and removed)
# PJSC20_infected_cell_barcodes_and_clusters_D1 <- subset(PJSC20_infected_cell_barcodes_and_clusters_updated, subset = PJSC20_infected_cell_barcodes_and_clusters_updated$Exp == 'Donor1' & PJSC20_infected_cell_barcodes_and_clusters_updated$percent.microsporidia > 2)
# cell.use <- PJSC20_infected_cell_barcodes_and_clusters_D1$...1
# Donor1.microsporidian <- subset(Donor1.microsporidian, cells = cell.use)
# PJSC20_clusters_D1 <- PJSC20_infected_cell_barcodes_and_clusters_D1$integrated_snn_res.0.2
# Donor1.microsporidian <- AddMetaData(Donor1.microsporidian, metadata = PJSC20_clusters_D1, col.name = "PJSC20_clusters")

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)

#Donor 2
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

#separate E. intestinalis transcripts
options(max.print = 100000)
rownames(gene.hashtag)
#Human transcripts = 1:36601
#Microsporidian transcripts = 36602:38612
Donor2.microsporidian <- subset(gene.hashtag[36602:38612,])
rownames(Donor2.microsporidian)
Donor2.microsporidian <- NormalizeData(Donor2.microsporidian, normalization.method = "LogNormalize", scale.factor = 10000)
Donor2.microsporidian <- FindVariableFeatures(Donor2.microsporidian, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Donor2.microsporidian)
Donor2.microsporidian <- ScaleData(Donor2.microsporidian, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
Donor2.microsporidian[["HTO"]] <- CreateAssayObject(counts = raw.HTO.count)
Donor2.microsporidian <- MULTIseqDemux(Donor2.microsporidian, assay = "HTO")
table(Donor2.microsporidian$MULTI_ID)
#12hpi     1dpi     3hpi     ctrl  Doublet Negative 
#3483     3219     3434      415      628     1100
remove(gene.hashtag)
Donor2.microsporidian <- subset(Donor2.microsporidian, idents = c("ctrl","3hpi","12hpi","1dpi"))
dim(Donor2.microsporidian) #10,551 cells

#### >>> FIX: subset infected cells directly (replace PJSC15/PJSC20 usage)
Donor2.microsporidian <- subset(Donor2.microsporidian, subset = percent.microsporidia > 2)
dim(Donor2.microsporidian)  # expect ~1k cells
#### <<<

#### (old code removed)
# cell.use <- PJSC15_infected_cell_barcodes_and_clusters_D2$...1
# Donor2.microsporidian <- subset(Donor2.microsporidian, cells = cell.use)
# PJSC20_clusters_D2 <- PJSC20_infected_cell_barcodes_and_clusters_D2$integrated_snn_res.0.2
# Donor2.microsporidian <- AddMetaData(Donor2.microsporidian, metadata = PJSC20_clusters_D2, col.name = "PJSC20_clusters")

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)

#Donor3 and Donor4
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

#separate E. intestinalis transcripts
options(max.print = 100000)
rownames(gene.hashtag)
#Human transcripts = 1:36601
#Microsporidian transcripts = 36602:38612
Donor3_4.microsporidian <- subset(gene.hashtag[36602:38612,])
rownames(Donor3_4.microsporidian)
Donor3_4.microsporidian <- NormalizeData(Donor3_4.microsporidian, normalization.method = "LogNormalize", scale.factor = 10000)
Donor3_4.microsporidian <- FindVariableFeatures(Donor3_4.microsporidian, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Donor3_4.microsporidian)
Donor3_4.microsporidian <- ScaleData(Donor3_4.microsporidian, features = all.genes)
raw.HTO.count <- NormalizeData(raw.HTO.count, normalization.method = "CLR")
Donor3_4.microsporidian[["HTO"]] <- CreateAssayObject(counts = raw.HTO.count)
Donor3_4.microsporidian <- MULTIseqDemux(Donor3_4.microsporidian, assay = "HTO")
table(Donor3_4.microsporidian$MULTI_ID)
#D3-1dpi  D3-2dpi  D3-ctrl  D4-1dpi  D4-2dpi  D4-ctrl  Doublet Negative 
#3500     5984      830     5583     4431     1572     6535     1920 

Donor3.microsporidian <- subset(Donor3_4.microsporidian, idents = c("D3-ctrl","D3-1dpi","D3-2dpi"))
Donor4.microsporidian <- subset(Donor3_4.microsporidian, idents = c("D4-ctrl","D4-1dpi","D4-2dpi"))

#### >>> FIX: subset infected cells directly for D3 and D4
Donor3.microsporidian <- subset(Donor3.microsporidian, subset = percent.microsporidia > 2)
dim(Donor3.microsporidian)  # ~1.9k cells expected

Donor4.microsporidian <- subset(Donor4.microsporidian, subset = percent.microsporidia > 2)
dim(Donor4.microsporidian)  # ~2.6k cells expected
#### <<<

#### (old code removed)
# cell.use <- PJSC15_infected_cell_barcodes_and_clusters_D3$...1
# Donor3.microsporidian <- subset(Donor3.microsporidian, cells = cell.use)
# PJSC20_clusters_D3 <- PJSC20_infected_cell_barcodes_and_clusters_D3$integrated_snn_res.0.2
# Donor3.microsporidian <- AddMetaData(Donor3.microsporidian, metadata = PJSC20_clusters_D3, col.name = "PJSC20_clusters")
# cell.use <- PJSC15_infected_cell_barcodes_and_clusters_D4$...1
# Donor4.microsporidian <- subset(Donor4.microsporidian, cells = cell.use)
# PJSC20_clusters_D4 <- PJSC20_infected_cell_barcodes_and_clusters_D4$integrated_snn_res.0.2
# Donor4.microsporidian <- AddMetaData(Donor4.microsporidian, metadata = PJSC20_clusters_D4, col.name = "PJSC20_clusters")

remove(raw.data)
remove(gene.hashtag)
remove(gene.singlet)
remove(percent.mt)
remove(raw.gene.count)
remove(raw.HTO.count)
remove(joint.bcs)
remove(all.genes)
remove(percent.microsporidia)
remove(Donor3_4.microsporidian)
# remove(PJSC20_clusters_D1)
# remove(PJSC20_clusters_D3)
# remove(PJSC20_clusters_D4)
remove(cell.use)

#Assign donor name into each cell and store the data in 'Exp' slot
Donor1.microsporidian$Exp <- "Donor1"
Donor2.microsporidian$Exp <- "Donor2"
Donor3.microsporidian$Exp <- "Donor3"
Donor4.microsporidian$Exp <- "Donor4"

saveRDS(Donor1.microsporidian, "Donor1.microsporidian.rds")
saveRDS(Donor2.microsporidian, "Donor2.microsporidian.rds")
saveRDS(Donor3.microsporidian, "Donor3.microsporidian.rds")
saveRDS(Donor4.microsporidian, "Donor4.microsporidian.rds")

#2022-05-20
#Perform integration and batch correction
#First, select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(Donor1.microsporidian, Donor2.microsporidian, Donor3.microsporidian, Donor4.microsporidian), nfeatures = 3000)
options(max.print = 100000)
features

#Identify anchors which will be used for integration
combined.anchors <- FindIntegrationAnchors(object.list = list(Donor1.microsporidian, Donor2.microsporidian, Donor3.microsporidian, Donor4.microsporidian), anchor.features = features)

#### >>> FIX: remove deletes of objects you never created
# remove(PJSC15_infected_cell_barcodes_and_clusters_D1)
# remove(PJSC15_infected_cell_barcodes_and_clusters_D2)
# remove(PJSC15_infected_cell_barcodes_and_clusters_D3)
# remove(PJSC15_infected_cell_barcodes_and_clusters_D4)
#### <<<

remove(combined.dataset)

#Combined datasets
combined.dataset <- IntegrateData(anchorset = combined.anchors)
dim(combined.dataset)
#6,908 cells with 1,990 features

#Next, let's see the integration results
DefaultAssay(combined.dataset) <- "integrated"
combined.dataset <- ScaleData(combined.dataset, verbose = FALSE)
combined.dataset <- RunPCA(combined.dataset)
##The following 3 features requested have zero variance (running reduction without them): Eint-010010, Eint-050010, Eint-070830
percent.pc <- combined.dataset@reductions$pca@stdev / sum(combined.dataset@reductions$pca@stdev) * 100
cum.percent.pc <- cumsum(percent.pc)
cum.percent.pc
#PC1-44

#Run UMAP
set.seed(2022)
combined.dataset <- RunUMAP(combined.dataset, dims = 1:44)
Idents(combined.dataset) <- "Exp"

png("PJSC23_UMAP_microsporidia_all_donors.png", width = 12, height = 8, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", pt.size = 0.1, cols = "Set1")
dev.off()

png("PJSC23_UMAP_microsporidia_all_donors_split.png", width = 40, height = 10, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", pt.size = 0.1, split.by = "Exp", ncol = 4, cols = "Set1")
dev.off()

#Let's see how %microsporidia are distributed
png("PJSC23_UMAP_percent.microsporidia_NoLabel.png", width = 12, height = 8, units = "cm", res = 300)
FeaturePlot(combined.dataset, features = "percent.microsporidia")
dev.off()

#Split according to the infection timepoints
combined.dataset$MULTI_ID <- recode_factor(combined.dataset$MULTI_ID, "D3-ctrl" = "ctrl", "D3-1dpi" = "1dpi", "D3-2dpi" = "2dpi", "D4-ctrl"  = "ctrl", "D4-1dpi" = "1dpi", "D4-2dpi" = "2dpi")
combined.dataset$MULTI_ID <- factor(combined.dataset$MULTI_ID, levels = c("ctrl","3hpi","12hpi","1dpi","2dpi","3dpi"))

png("PJSC23_UMAP_microsporidia_all_donors_split2.png", width = 32, height = 20, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "MULTI_ID", cols = "Set1", pt.size = 0.1) + facet_wrap(~combined.dataset$MULTI_ID + combined.dataset$Exp, ncol = 6)
dev.off()

png("PJSC16_UMAP_microsporidia_all_donors_split3.png", width = 24, height = 15, units = "cm", res = 300)
DimPlot(combined.dataset, reduction = "umap", split.by = "MULTI_ID", ncol = 3, cols = "Set1", pt.size = 0.1) 
dev.off()

#Perform the clustering
set.seed(2022)
Idents(combined.dataset) <- "orig.ident"
combined.dataset <- FindNeighbors(combined.dataset, dims = 1:44)

#Resolution 0.35 --> DECIDED TO USE 0.35
combined.dataset <- FindClusters(combined.dataset, resolution = 0.35)
png("PJSC23_UMAP_res0.35_No_label.png", width = 12, height = 8, units = "cm", res = 300)
DimPlot(combined.dataset, pt.size = 0.1, label = T, label.size = 4)
dev.off()

DimPlot(combined.dataset, split.by = "MULTI_ID", pt.size = 0.1) + facet_wrap(~combined.dataset$MULTI_ID + combined.dataset$Exp, ncol = 6)

#### >>> FIX: guard optional human cluster plots (only if metadata exists)
if ("PJSC20_clusters" %in% colnames(combined.dataset@meta.data)) {
  png("PJSC23_Human_clusters_mapped_res0.35.png", width = 12, height = 8, units = "cm", res = 300)
  DimPlot(combined.dataset, group.by = "PJSC20_clusters", pt.size = 0.4, label.size = 5)
  dev.off()
  
  png("PJSC23_human_clusters_split.png", width = 34, height = 15, units = "cm", res = 300)
  DimPlot(combined.dataset, group.by = "PJSC20_clusters", pt.size = 0.4, label.size = 5, split.by = "MULTI_ID") 
  dev.off()
}
#### <<<

png("PJSC23_UMAP_Split_Time_No_label_2.png", width = 25, height = 15, units = "cm", res = 300)
DimPlot(combined.dataset, split.by = "MULTI_ID", pt.size = 0.1, ncol = 3)
dev.off()

######################################
# (rest of your analysis continues unchanged)
