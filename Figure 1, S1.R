# HER2-positive project
# Load packages
library(Seurat)
library(DoubletFinder)
library(magrittr)
library(hdf5r)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
library(stringr)

# Load datasets
countpath <- "/home/chywang/Program/Singlecell_HER2/Data"
combine <- lapply(ID, function(i){CreateSeuratObject(counts = Read10X(data.dir = i),
                                                     min.cells = 3, 
                                                     min.features = 800,
                                                     project = str_split_fixed(i, "/", 9)[ ,8])})
combine <- lapply(combine, function(i){subset(i, subset = `nFeature_RNA` <= 6000)})
combine <- lapply(combine, function(i){
  PercentageFeatureSet(i, pattern = "^MT-", col.name = "percent.mt")})
combine <- lapply(combine, function(i){subset(i, percent.mt < 10)})

seuratdeal <- function(SeuratObject, resolution){
  SeuratObject <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
  # Find Variable Features
  SeuratObject <- FindVariableFeatures(SeuratObject, selection.method = "vst", nfeatures = 2000)
  SeuratObject <- ScaleData(SeuratObject)
  SeuratObject <- RunPCA(SeuratObject, verbose = F, dims = 1:30, npcs = 30)
  SeuratObject <- RunUMAP(SeuratObject, dims=1:16)
  SeuratObject <- FindNeighbors(SeuratObject, dims = 1:16) %>% FindClusters(resolution = resolution)
}
combine <- lapply(combine, function(i){seuratdeal(i, 0.7)})
doubletratio <- c()
for (i in 1:length(combine)) {
  doubletratio[i] <- (combine[[i]]@meta.data %>% nrow()) / 500 * 0.004}

# Delete doublet ----
DelDoublet <- function(SeuratObject, DoubletRate){
  sweep.res.list <- paramSweep_v3(SeuratObject, PCs = 1:16, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  homotypic.prop <- modelHomotypic(SeuratObject$seurat_clusters)   # celltype is the best choice
  nExp_poi <- round(DoubletRate*ncol(SeuratObject)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  SeuratObject <- doubletFinder_v3(SeuratObject, PCs = 1:16, pN = 0.25, pK = pK_bcmvn, 
                                   nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
  return(SeuratObject)
}
combine <- mapply(DelDoublet, combine, doubletratio)
# parallel
for (i in 1:length(combine)) {
  names(combine[[i]]@meta.data)[combine[[i]]@meta.data %>% ncol()] <- "Doubletstate"
}
combine <- lapply(combine, function(i){subset(i, Doubletstate == "Singlet")})
saveRDS(combine, "Combine.rds")

# Lineages markers
CombineMarkers <- FindAllMarkers(combine, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
DotPlot(combine, features = c("CD2", "CD3D",
                              "FGFBP2", "FCGR3A",
                              "MS4A1", "CD19",
                              "SDC1", "TNFRSF17",
                              "CD14", "S100A9",
                              "TPSAB1", "CPA3"), group.by = "seurat_clusters")

