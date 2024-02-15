
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggraph)
library(limma)
library(magrittr)
library(tibble)
library(tidyverse)
library(SeuratObject)
library(clustree)
library(SingleR)
library(celldex)
#
setwd("~/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/11.MMsub_copykat")
load("~/cellranger.7.0.1/2.MM/500_5000_10/11.copykat/MM.Rda")
load("~/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/9.non_mm/Melanoma.name.Rda")
table(Idents(MM))
table(MM$orig.ident)
DefaultAssay(object = MM) <- "integrated"
scRNA <- subset(MM, cells = Melanoma.name)
scRNA <- subset(scRNA, idents =  "aneuploid")
cellfile <- "MMsub"
scRNA <- RunPCA(scRNA, verbose = FALSE)
scRNA <- RunUMAP(scRNA, dims = 1:50, verbose = FALSE)
scRNA <- FindNeighbors(scRNA, dims = 1:50, verbose = FALSE)
scRNA <- FindClusters(scRNA, verbose = FALSE, resolution = 0.1)
save(scRNA, file = paste0(cellfile, "_dim50_r01.Rda"))

Idents(scRNA) <- scRNA$integrated_snn_res.0.1
P1 <- DimPlot(scRNA, label = TRUE) + NoLegend()
pdf(paste0("P1_", cellfile, "_umap_dim50_r01.pdf"), width = 8, height = 6)
P1
dev.off()

dir.create("3.Normalization")
DefaultAssay(object = scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)         
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), features = all.genes) 
save(scRNA, file = paste0("3.Normalization/", cellfile, "_exmt_cc_NormalizeData.Rda"))



