
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
library(copykat)

setwd("~/cellranger.7.0.1/2.MM/500_5000_10/11.copykat")

load("~/cellranger.7.0.1/2.MM/500_5000_10/8.celltype/MM_celltype.Rda")
DefaultAssay(MM) <- "RNA"
table(Idents(MM))

split_seurat <- SplitObject(MM, split.by = "orig.ident")
samples <- c("MM2", "MM6", "MM7", "MM13", "MM14", "MM21", "MM24", "MM25", 
             "MM27", "MM28", "MM29", "MM30", "MM33", "MM35", "MM38", "MM40", "MM47", "MM48", "MM49", "MM50")
copykat.test.list <- list()
for (i in samples) {
  
  exp.rawdata <-GetAssayData(split_seurat[[i]], assay = "RNA", slot = "count")
  copykat.test <- copykat(rawmat=exp.rawdata, 
                          id.type="S", 
                          ngene.chr=5, 
                          win.size=25, 
                          KS.cut=0.1, 
                          sam.name="test", 
                          distance="euclidean", 
                          n.cores=16)
  save(copykat.test, file = paste0(i,"copykat.test.Rda"))
  copykat.test.list[[i]] <- copykat.test
  
}
save(copykat.test.list, file = "copykat.test.list.Rda")
