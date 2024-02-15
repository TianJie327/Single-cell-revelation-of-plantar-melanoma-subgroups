
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
setwd("~/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/11.MMsub_copykat/20.MapQuery")
#
load("~/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/11.MMsub_copykat/14.melanoma_Diploid_melanocytes_schwann/Melanocytes.Rda")
load("~/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/11.MMsub_copykat/14.melanoma_Diploid_melanocytes_schwann/Schwann.Rda")
load("~/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/11.MMsub_copykat/12.GSEA_forloop/MM_celltype.Rda")


MM.query <-merge( x = melanocytes, y = c(schwann), merge.data = TRUE)
MM.integrated <- MM

MM.integrated <- RunUMAP(MM, dims = 1:50, reduction = "pca", return.model = TRUE)
p1=DimPlot(MM.integrated,label = T)
p1
pdf("MM.integrated.pdf", width = 8, height = 6)
p1
dev.off()



MM.anchors <- FindTransferAnchors(reference = MM.integrated, query = MM.query,
                                  dims = 1:50, 
                                  reference.reduction = "pca")

MM.query <- TransferData(anchorset = MM.anchors, reference = MM.integrated, query = MM.query,
                         refdata = list(celltype = "cellType"))
MM.query <- IntegrateEmbeddings(anchorset = MM.anchors, reference = MM.integrated,
                                query = MM.query, new.reduction.name = "ref.pca")
MM.query <- ProjectUMAP(query = MM.query, query.reduction = "ref.pca", reference = MM.integrated,
                        reference.reduction = "pca",
                        reduction.model = "umap"
) 


p1 <- DimPlot(MM.integrated, reduction = "pca", group.by = "cellType", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(MM.query, reduction = "ref.pca", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf("ref.pca.pdf", width = 10, height = 6)
p1 + p2
dev.off()

p1 <- DimPlot(MM.integrated, reduction = "umap", group.by = "cellType", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(MM.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
pdf("ref.umap.pdf", width = 12, height = 6)
p1 + p2
dev.off()


