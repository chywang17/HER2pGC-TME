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

# Filter B cells
Bcell <- subset(Combine, CellType == "B cell")

# Run the standard workflow for visualization and clustering
Bcell <- NormalizeData(Bcell, normalization.method = "LogNormalize", scale.factor = 1e4) 
Bcell <- FindVariableFeatures(Bcell, selection.method = "vst", nfeatures = 2000)
Bcell <- ScaleData(Bcell)
Bcell <- RunPCA(Bcell, verbose = FALSE)
Bcell <- RunHarmony(Bcell, "orig.ident", plot_convergence = TRUE, lambda = 1.5)
Bcell <- RunUMAP(Bcell, reduction = "harmony", dims = 1:30)
Bcell <- FindNeighbors(Bcell, reduction = "harmony", dims = 1:30)
Bcell <- FindClusters(Bcell, resolution = 1)
DimPlot(Bcell, group.by = "seurat_clusters", label = T)

# Ro/e analysis
data <- RoeIndex(Bcell, "Group", "Celltype")
data <- data[ ,-3:-4] # delete Freq, Expected 
data <- reshape2::acast(data, Var1 ~ Var2) %>% t()
bk <- c(seq(0.5, 1, by = 0.01))
Group <- data.frame(Stage = c("Pre-treatment", "Post-treatment", "Pre-treatment", "Post-treatment"),
                    Group = c("Control", "Control", "HER2", "HER2"))
rownames(Group) <- colnames(data)
Group$Group <- factor(Group$Group)
Group$Stage <- factor(Group$Stage)

ann_colors = list(
  Group = c(HER2 = "#D3BA68FF", Control = "#65A479FF"),
  Stage = c(`Pre-treatment` = "#D5695DFF", `Post-treatment` = "#5D8CA8FF"))
pl <- pheatmap(data, scale="none",
               cluster_rows = T,
               cluster_cols = F,
               show_rownames = T,
               show_colnames = F,
               color = paletteer::paletteer_c("grDevices::Purple-Orange", n = length(bk), direction = -1),
               annotation_colors = ann_colors,
               border=NA,
               fontsize = 10,
               breaks = bk,
               display_numbers = TRUE,
               annotation_col = Group,
               number_color = "white",
               fontsize_number = 10,
               legend_breaks = seq(0.5, 1, 0.5),
               gaps_col = 2)


# Load datasets (TPM and Clinical datasets)
HER2TPM <- readRDS("~/Data/Program/sc_HER2GC/2_Code/RCode/HER2TPM.rds")
Clinical <- read.table("./Data/SYSUCC_Bulk.txt", header = T)
Clinical$SampleID <- paste0("A", Clinical$SampleID)
Clinical <- column_to_rownames(Clinical, "SampleID")

# Run GSVA
enrichscore <- gsva(as.matrix(HER2TPM), list(ERBB2 = ERBB2,
                                             Bcell = Bcell,
                                             Tfh = Tfh,
                                             Tex = Tex,
                                             Exh = Exh),
                    method= "ssgsea") %>% as.data.frame()
enrichscore <- t(enrichscore) %>% as.data.frame()
# Combine with clinical data
enrichscore <- merge(enrichscore, Clinical, by = "row.names")
# Correlation plot
ggplot(enrichscore, aes(x = ERBB2, y = Exh))+ 
  geom_smooth(method = "lm")+
  ggpubr::stat_cor(method = "spearman")+
  geom_point(size = 0.8)+
  theme_bw()

# TCGA-STAD Cohort
TCGA_STAD <- readRDS("~/Data/Program/sc_HER2GC/2_Code/RCode/Data/TCGA_STAD_TPM.rds")
TCGA_STAD <- TCGA_STAD[ ,-which(str_split_fixed(colnames(TCGA_STAD), "[.]", 4)[ ,4] == "11A")]
enrichscore <- gsva(as.matrix(TCGA_STAD), list(ERBB2 = ERBB2,
                                               Bcell = Bcell,
                                               Tfh = Tfh,
                                               Tex = Tex,
                                               Exh = Exh),
                    method= "ssgsea") %>% as.data.frame()
enrichscore <- t(enrichscore) %>% as.data.frame()
enrichscore$Group <- "HER2+"
enrichscore$Group[enrichscore$ERBB2 <= median(enrichscore$ERBB2)] <- "HER2-"

# Group
ColorPanel <- c("#93C6E1FF", "#F8B150FF")
names(ColorPanel) <- c("HER2-", "HER2+")
ggplot(enrichscore, aes(x = Group, y = Tfh, fill = Group))+
  stat_boxplot(geom = "errorbar",position=position_dodge(width = 0.2), width = 0.1)+
  geom_boxplot(aes(fill = Group), position=position_dodge(width = 0.2), width = 0.4)+
  stat_compare_means(comparisons = list(c("HER2+", "HER2-")))+
  ylab("Tfh")+
  xlab("Group")+
  scale_fill_manual(values = ColorPanel)+
  ggprism::theme_prism()

# Proteomics
HER2Pro <- read.delim("/Users/chaoyewang/Data/Program/sc_HER2GC/2_Code/RCode/Data/FDR.txt", header = T)
HER2Pro <- HER2Pro[ ,c(-1,-3,-4)]
HER2Pro <- aggregate(HER2Pro[ ,-1], list(HER2Pro[ ,1]), mean)
HER2Pro <- column_to_rownames(HER2Pro, "Group.1")
enrichscore <- gsva(as.matrix(HER2Pro), 
                    list(ERBB2 = ERBB2,
                         Bcell = Bcell,
                         Tfh = Tfh,
                         Tex = Tex),
                    method= "ssgsea") %>% 
  t() %>% 
  as.data.frame()
enrichscore$PDL1 <- as.numeric(HER2Pro["CD274", ])
enrichscore$ERBB2Gene <- as.numeric(HER2Pro["ERBB2", ])
ClinicalData <- read.delim("/Users/chaoyewang/Data/Program/sc_HER2GC/2_Code/RCode/Data/Clinical.txt", header = T)
ClinicalData <- column_to_rownames(ClinicalData, "Firmiana.ID")
enrichscore <- merge(enrichscore, ClinicalData, by = "row.names")
enrichscore <- subset(enrichscore, !is.na(HER2.status))

plGroupColor <- c("#93C6E1FF", "#F8B150FF")
names(plGroupColor) <- c("Negative", "Positive")

ggplot(enrichscore, aes(x = HER2.status, y = Bcell, fill = HER2.status))+
  stat_boxplot(geom = "errorbar",position=position_dodge(width = 0.2), width = 0.1)+
  geom_boxplot(aes(fill = HER2.status), position=position_dodge(width = 0.2), 
               width = 0.4)+
  stat_compare_means(comparisons = list(c("Negative", "Positive")))+
  scale_fill_manual(values = plGroupColor)+
  ggprism::theme_prism()

CD27Ratio <- read.table("./Data/mIHCData/CD27Ratio.txt", header = T)
CD27Ratio$Stage <- str_split_fixed(CD27Ratio$Group, "-", 2)[ ,2]

ColorPalette <- c("#3D5DA9", "#B0395D")
names(ColorPalette) <- c("pre", "pro")

ggplot(CD27Ratio, aes(x = Stage, y = CD11CB, color = Stage))+
  stat_boxplot(geom = "errorbar",position=position_dodge(width = 0.4), width = 0.2)+
  geom_boxplot(aes(color = Stage), position=position_dodge(width = 0.4), 
               width = 0.6, outlier.shape = NA)+
  geom_jitter(aes(color = Stage), width = 0.1, alpha = 0.4)+
  stat_compare_means(comparisons = list(c("pre", 
                                          "pro")), method = "t")+
  scale_color_manual(values = ColorPalette)+
  ylab("CD11c+ B cells")+
  xlab("Stage")+
  ggprism::theme_prism()+
  theme(legend.position = 'none')
aov(CD27B ~ Stage,
    data = CD27Ratio) %>% 
  summary()


CD27Dis <- read.table("./Data/mIHCData/CD27Dis.txt", header = T)
CD27Dis$Stage <- str_split_fixed(CD27Dis$Group, "-", 2)[ ,2]

ggplot(CD27Dis, aes(x = CD11CBDis, fill = Stage)) + 
  geom_density(alpha = 0.5)+
  theme(legend.position = 'none')+
  labs(title = "Density Distribution", x = "Value", y = "Density")

Pre_CD27Dis <- subset(CD27Dis, Stage == "pro")
Pre_CD27Dis_Frame <- data.frame(Distance = c(Pre_CD27Dis$CD27BDis, Pre_CD27Dis$CD11CBDis))
Pre_CD27Dis_Frame$Celltype <- "CD11C"
Pre_CD27Dis_Frame$Celltype[1:43] <- "CD27"

# Distance to ERBB2+ cancer cells
ColorPalette <- c("#445B85", "#CB7767")
names(ColorPalette) <- c("CD11C", "CD27")

ggplot(Pre_CD27Dis_Frame, aes(x = Celltype, y = Distance, color = Celltype))+
  stat_boxplot(geom = "errorbar",position=position_dodge(width = 0.4), width = 0.2)+
  geom_boxplot(aes(color = Celltype), position=position_dodge(width = 0.4), 
               width = 0.6, outlier.shape = NA)+
  geom_jitter(aes(color = Celltype), width = 0.1, alpha = 0.4)+
  stat_compare_means(comparisons = list(c("CD11C", "CD27")), method = "t")+
  scale_color_manual(values = ColorPalette)+
  ylab("Distance (µm) to HER2(+) cancer cells")+
  xlab("Cell types")+
  ggprism::theme_prism()+
  theme(legend.position = 'none')