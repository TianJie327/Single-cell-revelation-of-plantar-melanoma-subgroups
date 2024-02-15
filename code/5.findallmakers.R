
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggraph)
library(magrittr)
library(tibble)
library(tidyverse)
library(SeuratObject)
load("~/cellranger.7.0.1/2.MM/500_5000_10/20.article_figs/4.fig1/MM_all.celltype.article.Rda")
setwd("~/cellranger.7.0.1/2.MM/500_5000_10/20.article_figs/11.supplementary/table/2.DEG/fig1")
MM <- MM_all
table(Idents(MM))
DefaultAssay(MM) <- "RNA"
markers <- FindAllMarkers(MM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
write.csv(markers, file = "Allmarker.csv")
save(markers, file = "Allmarker.Rda")
