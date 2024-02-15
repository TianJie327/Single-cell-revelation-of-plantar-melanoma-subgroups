
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(NMF)

load("/gdata01/user/tianjie/cellranger.7.0.1/2.MM/500_5000_10/20.article_figs/4.fig1/MM_adjust_CAFs_20231206.Rda")

DefaultAssay(object = MM) <- "RNA"
Demo <- MM
cellchat <- createCellChat(Demo@assays$RNA@data, meta = Demo@meta.data,group.by = "celltype5")
cellchat <- setIdent(cellchat, ident.use = "celltype5") 
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use 

cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
save(cellchat, file = "cellchat.Rda")

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
save(cellchat, file = "cellchat.Rda")

df.net <- subsetCommunication(cellchat) 
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
save(cellchat, file = "cellchat.Rda")

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pdf("cellchat.pdf", height = 20, width = 20)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], vertex.label.cex = 0.5)
}
dev.off()

mat <- cellchat@net$weight
par(mfrow = c(1,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

netVisual_circle(cellchat@net$weight  , vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
gg1 <- netAnalysis_signalingRole_scatter(cellchat)


ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 15, height = 30, font.size = 8)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 15, height = 30, font.size = 8)
ht1 + ht2


save(cellchat, file = "cellchat.Rda")


