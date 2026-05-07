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

# CSFE
CSFE_NB_Pan_T <- read.table("./Data/mIHCData/CSFE_NB_Pan_T.txt", header = T)
CSFE_NB_Pan_T$Group <- factor(CSFE_NB_Pan_T$Group, levels = c("D3-BT", "D5-BT", "D3-BT+", "D5-BT+"))
ggplot(CSFE_NB_Pan_T, aes(Group, D5.CFSE))+
  stat_summary(mapping = aes(fill = Group),
               fun = mean,
               geom = "bar",
               fun.args = list(mult = 1), 
               width = 0.7)+
  stat_summary(fun.data = mean_se, 
               fun.args = list(mult = 1), 
               geom = "errorbar", 
               width=0.2)+
  geom_jitter(aes(fill = Group),
              position = position_jitter(0.1), 
              shape=21, 
              size = 2, 
              alpha=0.9)+
  stat_compare_means(comparisons = list(c("D3-BT", "D3-BT+"),
                                        c("D5-BT", "D5-BT+")),
                     method = "t.test")+
  scale_fill_manual(values = ColorPalette)+
  theme_bw()+
  theme(legend.position = "none")

CSFE_NB_Pan_T <- read.table("./Data/mIHCData/CSFE_NB_Pan_T_D5.txt", header = T)
ggplot(CSFE_NB_Pan_T, aes(Group, D5.CFSE))+
  stat_summary(mapping = aes(fill = Group),
               fun = mean,
               geom = "bar",
               fun.args = list(mult = 1), 
               width = 0.7)+
  stat_summary(fun.data = mean_se, 
               fun.args = list(mult = 1), 
               geom = "errorbar", 
               width = 0.2)+
  geom_jitter(aes(fill = Group),
              position = position_jitter(0.1), 
              shape=21, 
              size = 2, 
              alpha = 0.9)+
  stat_compare_means(comparisons = list(c("T", "T+"),
                                        c("T", "BT+"),
                                        c("T+", "BT+")),
                     method = "t.test")+
  scale_fill_manual(values = ColorPalette)+
  theme_bw()+
  theme(legend.position = "none")






