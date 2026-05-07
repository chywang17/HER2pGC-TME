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

# T cell
Tcell$treatmentstate <- paste(Tcell$treatment, ".", Tcell$state, sep = "")
Tcell <- subset(combinelist, celltype == "T cell")
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
# Find Variable Features
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- ScaleData(Tcell)
Tcell <- RunHarmony(Tcell, group.by.vars = "orig.ident", lambda = 2)
Tcell <- RunUMAP(Tcell, reduction = "harmony", dims = 1:30)
Tcell <- FindNeighbors(Tcell, reduction = "harmony", dims = 1:30) %>% FindClusters()
TcellMarkers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
saveRDS(Tcell, "Tcell_afterbatch,rds")
DimPlot(Tcell, group.by = "seurat_clusters", label = T)+
  Tcellcolor(1)

# 3 Cell fraction analysis ----
# 3.1 Roe analysis ----
data <- RoeIndex(Tcell, "Status", "seurat_clusters")
data <- data[ ,-3:-4]
data$Var1 <- factor(data$Var1, levels = c("Pre.HER2(+)", "Pre.HER2(-)", "Post.HER2(+)", "Post.HER2(-)"))
data <- reshape2::acast(data, Var1 ~ Var2) %>% t()
bk <- c(seq(0, 1.5, by = 0.1))
pl <- pheatmap(data, scale="none",
               cluster_rows = F,
               cluster_cols = F,
               show_rownames = T,
               show_colnames = T,
               color = rev(paletteer_c("grDevices::Purple-Orange", n = length(bk))),
               border=NA,
               fontsize = 10,
               breaks = bk,
               display_numbers = TRUE,
               #annotation_col = group,
               number_color = "white",
               fontsize_number = 10,
               legend_breaks = seq(0, 1, 0.1))

# STARTRAC analysis
CD8_Metadata <- readRDS("~/Data/Program/sc_HER2GC/2_Code/RCode/Data/CD8_Metadata.rds")
CD4_Metadata <- readRDS("~/Data/Program/sc_HER2GC/2_Code/RCode/Data/CD4_Metadata.rds")
PostID <- c("aPA7", "aPA8", "aPA11", "aPA9", "aPA10", "aPA12")
CD8_Metadata <- subset(CD8_Metadata, Frequency > 1 & orig.ident %in% PostID)
CD4_Metadata <- subset(CD4_Metadata, Frequency > 1 & orig.ident %in% PostID)

# Prepare input files
# CD8 files
CD8_TCR_Input <- data.frame(Cell_Name = rownames(CD8_Metadata))
CD8_TCR_Input <- cbind(CD8_TCR_Input,
                       CD8_Metadata[ ,c("CTaa", "orig.ident", "Patient", 
                                        "Group", "orig.ident", "Celltype", 
                                        "orig.ident")])
colnames(CD8_TCR_Input) <- c("Cell_Name", "clone.id", "clone.status", 
                             "patient", "sampleType", "stype", "majorCluster", 
                             "loc")
CD8_TCR_Input$clone.status <- "Clonal"
CD8_TCR_Input$stype <- "CD8"
CD8_TCR_Input$loc <- "T"

# Run Startrac
Out <- Startrac.run(CD8_TCR_Input, proj = "HER2", verbose = T)
ExpanOut <- subset(Out@cluster.data, aid != "HER2")
ExpanOut$CliType <- "HER2p"
ExpanOut$CliType[ExpanOut$aid %in% c("PA7", "PA8", "PA11")] <- "HER2n"

pldata <- subset(ExpanOut, majorCluster == "CD8.c11.Tex(CXCL13)")
ggplot(pldata, aes(x = Group, y = expa, fill = Group))+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 1), width = 0.4)+
  geom_boxplot(aes(fill = Group), position = position_dodge(width = 1), 
               width = 1)+
  stat_compare_means(comparisons = list(c("neg_Pre", "neg_Post"), 
                                        c("pos_Pre", "pos_Post")))+
  ylab("Expansion score")+
  xlab("Group")+
  ggprism::theme_prism()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
