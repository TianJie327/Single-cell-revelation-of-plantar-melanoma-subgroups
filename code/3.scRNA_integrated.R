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

setwd("~/cellranger.7.0.1/2.MM")

dir.create("1.inter")
  
samples <- c("MM2", "MM6", "MM7", "MM13", "MM14", "MM21", "MM24", "MM25", 
             "MM27", "MM28", "MM29", "MM30", "MM33", "MM35", "MM38", "MM40", "MM47", "MM48", "MM49", "MM50"
             )

for (i in samples){
        seurat_data <- Read10X(data.dir = paste0("/gdata01/user/tianjie/cellranger.7.0.1/soupx/", i, "/", i, "/outs/", i ,"_strainedCounts_0.2/"))
        doublets <- read.table(paste0("/gdata01/user/tianjie/cellranger.7.0.1/20221016/", i, "/", i, "/outs/filtered_feature_bc_matrix/", i, "_doublet.txt"), header = T, sep = ",", row.names = 1)
        colnames(doublets) <- c("Doublet_score", "Is_doublet")
        seurat_obj <- CreateSeuratObject(counts = seurat_data,    
                                         meta.data = doublets,
                                         min.cells = 3, 
                                         min.features = 200,
                                         project = i)
        seurat_obj@meta.data$barcode <- row.names(seurat_obj@meta.data)
        seurat_obj <- RenameCells(seurat_obj, new.names = paste0(i, "_", 1:ncol(seurat_obj)))                                 
        assign(i, seurat_obj)
}

MM <- merge(x = MM2, y = c(MM6, MM7, MM13, MM14, MM21, MM24, MM25, MM27, MM28, MM29, MM30, MM33, MM35, MM38, MM40, MM47, MM48, MM49, MM50
                           ), merge.data = TRUE)

table(MM$Is_doublet)

MM <- subset(MM, subset = Is_doublet == "False")
table(MM$Is_doublet)


table(grepl("^MT-", row.names(MM)))
MM <- PercentageFeatureSet(MM, pattern = "^MT-", col.name = "percent.mt")
MM <- subset(MM, subset = nFeature_RNA >= 500 & nFeature_RNA <= 5000 & percent.mt <= 10)

HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
HB_m <- match(HB.genes, rownames(MM@assays$RNA))
HB.genes <- rownames(MM@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
MM[["percent.HB"]] <- PercentageFeatureSet(MM, features = HB.genes)

col.num <- length(levels(MM@active.ident))
violin <- VlnPlot(MM,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"),
        cols     = rainbow(col.num),
        pt.size  = 0.01,
        ncol     = 1
) +
        theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
        )


pdf("1.inter/P1_QC_VlnPlot.pdf", width = 10, height = 32) 
violin
dev.off()

p1 <- FeatureScatter(MM, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(MM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p3 <- FeatureScatter(MM, feature1 = "nCount_RNA", feature2 = "percent.HB")
P2 <- p1 | p2 | p3

pdf("1.inter/P2_QC_FeatureScatter.pdf", width = 20)
P2
dev.off()
save(MM, file = "1.inter/1_MM_before_integ.Rda")


split_seurat <- SplitObject(MM, split.by = "orig.ident")
rm(MM)
split_seurat <- split_seurat[samples]
for (i in 1:length(split_seurat)) {
        split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
        split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
        split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), conserve.memory = TRUE)
}
#save(split_seurat, file = "1.inter/3_split_seurat.Rda")
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)
#save(integ_features, file = "1.inter/3_1_sinteg_features.Rda")
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)
#save(split_seurat, file = "1.inter/3_2_split_seurat.Rda")
split_seurat <- lapply(X = split_seurat, FUN = RunPCA, features = integ_features)

integ_anchors <- FindIntegrationAnchors(
        object.list = split_seurat,
        reference = c(1, 2), #MM2, MM6 
        normalization.method = "SCT",
        anchor.features = integ_features,
        dims = 1:30,
        reduction = "rpca",
        k.anchor = 20
)
rm(split_seurat)
#save(integ_anchors, file = "1.inter/3_3_integ_anchors.Rda")
seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT", dims = 1:30)
save(seurat_integrated, file = "1.inter/4_integrated_seurat.Rda")


seurat_integrated <- RunPCA(object = seurat_integrated)
P3 <- PCAPlot(seurat_integrated, split.by = "orig.ident", ncol = 1)
#P3
pdf("1.inter/P3_PCAPlot_seurat_integrated.pdf", width = 10, height = 8*length(samples))
P3
dev.off()

seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:50, reduction = "pca")
P4 <- DimPlot(seurat_integrated)
#P4
pdf("1.inter/P4_umap.pdf", width = 10)
P4
dev.off()

P5 <- DimPlot(seurat_integrated, split.by = "orig.ident", ncol = 1)
#P5
pdf("1.inter/P5_umap_split.pdf", width = 10, height = 8*length(samples))
P5
dev.off()

MM <- seurat_integrated
rm(seurat_integrated)
MM <- FindNeighbors(MM, dims = 1:50, reduction = "pca", verbose = FALSE)
MM <- FindClusters(MM, verbose = FALSE, resolution = 0.1)
P6 <- DimPlot(MM, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5)
P6
pdf("1.inter/P6_umap_MM_r0.1.pdf", width = 10)
P6
dev.off()

save(MM, file = "1.inter/5_MM.Rda")

