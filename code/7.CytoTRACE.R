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
library(CytoTRACE)

load("~/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/11.MMsub_copykat/4.clustree/MMsub.Rda")

data.name = 'Mel'
context_str = paste("CytoTRACE-", data.name, "-all_mode-", format(Sys.time(), "%b%d_%H_%M_%S"), "-", sep="") 
context_str

sce <- scRNA 
DefaultAssay(sce) <- "RNA"
sce$cellType <- Idents(sce)
table(sce$cellType)
phe <- sce$cellType
phe = as.character(phe)
names(phe) <- rownames(sce@meta.data)
mat <- as.matrix(sce@assays$RNA@counts)
results <- CytoTRACE(mat = mat,enableFast = T)
plotCytoGenes(results, numOfGenes = 10,outputDir = context_str)
emb <- sce@reductions$umap@cell.embeddings
plotCytoTRACE(results, phenotype = phe, gene = "MITF", emb = emb, outputDir = paste0('MITF-', context_str))


