# Load packages
library(ggplot2)
library(tidyverse)

# Load datasets
PDL1Stratify_Pan <- read.table("./PDL1Stratify_Pan.txt", header = T)
ColorPalette <- c("#0076C0FF", "#95C65CFF")
names(ColorPalette) <- c("PDL1Neg", "PDL1Pos")
ggplot(PDL1Stratify_Pan, aes(x = Group, y = Pan, color = Group))+
  stat_boxplot(geom = "errorbar",position=position_dodge(width = 0.2), width = 0.1)+
  geom_boxplot(aes(color = Group), position=position_dodge(width = 0.2), 
               width = 0.4)+
  geom_jitter(aes(color = Group), shape = 21, width = 0.1)+
  stat_compare_means(comparisons = list(c("PDL1Neg", 
                                          "PDL1Pos")))+
  scale_color_manual(values = ColorPalette)+
  ylab("Pantothenic acid level")+
  xlab("Group")+
  ggprism::theme_prism()
one.way <- aov(Pan ~ Group, data = PDL1Stratify_Pan)
summary(one.way)
