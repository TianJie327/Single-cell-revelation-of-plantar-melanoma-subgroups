################################################################################
library("Capybara")
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
setwd("~/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/11.MMsub_copykat/21.Capybara")

#########################################################################################

load("MM_celltype.Rda")

stone.et.al.final <- read.table("df1.csv", header = T, sep = ",", row.names = 1, stringsAsFactors = F)
# MM <- subset(x = MM, downsample = 500) 
# stone.et.al.final <- GetAssayData(object = MM, assay = "RNA", slot = "counts")
# stone.et.al.final <- as.data.frame(stone.et.al.final)
# stone.et.al.final <- round(stone.et.al.final)
# MM <- subset(x = MM, downsample = 500) 

library(stringr)

names <- rownames(stone.et.al.final)
rownames(stone.et.al.final) <- str_to_title(names)


#########################################################################################
# File path
bulk.raw.path <- system.file("extdata", "Bulk Reference Raw.Rds", package = "Capybara")
bulk.rpkm.path <- system.file("extdata", "Bulk Reference RPKM.Rds", package = "Capybara")
# Read the matrices
bulk.raw <- readRDS(bulk.raw.path)
head(bulk.raw)
bulk.rpkm <- readRDS(bulk.rpkm.path)

single.round.QP.analysis(bulk.raw, stone.et.al.final, scale.bulk.sc = "scale", unix.par = F, 
                         force.eq = 1, n.cores = 1, save.to.path = "./", 
                         save.to.filename = "stone_bulk_classification_qp")

save(stone.et.al.final, file = "stone.et.al.final.Rda")
#####################################################################


## Load QP results
qp.rslt <- read.csv("./stone_bulk_classification_qp_scale.csv", row.names = 1, header = T, stringsAsFactors = F)
head(qp.rslt)
## Reshape the data
qp.rslt.sub <- qp.rslt[,c(1:(ncol(qp.rslt) - 2))]

## Background matrix
background.qp.fpath <- system.file("extdata", "MCA Embryonic Background.Rds", package = "Capybara")
background.mca <- readRDS(background.qp.fpath)
background.mtx <- background.mca[[2]]

## Correlation Analysis
mtx.test <- t(qp.rslt.sub[, colnames(background.mtx)])
ref.test <- t(background.mtx)

## Pearson's Correlation Calculation
corr.mtx <- WGCNA::cor(ref.test, mtx.test)

## Setup a correlation cutoff to the 90th quantile of the correlation matrix
correlation.cutoff <- quantile(corr.mtx, 0.90)

## Binarization based on the correlation
new.corr.bin <- corr.mtx
new.corr.bin[which(new.corr.bin >= correlation.cutoff)] <- 1
new.corr.bin[which(new.corr.bin < correlation.cutoff)] <- 0
new.corr.bin <- as.data.frame(new.corr.bin)
save(new.corr.bin, file = "new.corr.bin.Rda")


###########################################################################################################

# Count
count.in.cat <- c()
unique.cat <- unique(unlist(lapply(strsplit(rownames(new.corr.bin), "_"), function(x) x[1])))
for (uc in unique.cat) {
  curr.subset <- new.corr.bin[which(startsWith(rownames(new.corr.bin), uc)), c(1:ncol(stone.et.al.final))]
  count.in.cat[uc] <- sum(colSums(curr.subset) >= nrow(curr.subset) * 0.80)
}

count.in.cat <- as.data.frame(count.in.cat)
count.in.cat$perc <- round(count.in.cat$count.in.cat *100/sum(count.in.cat$count.in.cat), digits = 3)
save(count.in.cat, file = "count.in.cat.Rda")
final.cell.types.fetal <- rownames(count.in.cat)[which(count.in.cat$count.in.cat > 100)]
save(final.cell.types.fetal, file = "final.cell.types.fetal.Rda")

############################################################################################################################################


# Background cells
mca <- read.csv("~/code_learn/2.Capybara/2.exmple/Capybara-master/inst/extdata/MCA_CellAssignments.csv",
                row.names = 1, header = T, stringsAsFactors = F)

# Read the meta data


mca.meta <- data.frame(row.names = mca$Cell.name, 
                       tissue = mca$Tissue,
                       cell.bc.tissue = unlist(lapply(strsplit(mca$Cell.name, "_"), function(x) x[1])),
                       cell.type = mca$Annotation,
                       stringsAsFactors = F)

cardiac.rp.all.meta <- mca.meta[which(mca.meta$cell.bc.tissue %in% final.cell.types.fetal), ]

mca.counts.all.involved <- NULL
tissues.to.read <- unique(cardiac.rp.all.meta$tissue)
tissues.to.read <- gsub("_", "", tissues.to.read) #文件夹里面没有这些符号
tissues.to.read <- gsub("-", "", tissues.to.read)
general.path <- "~/code_learn/2.Capybara/5.Capybara/MCA_BatchRemove_dge/rmbatch_dge/"
for (i in 1:length(tissues.to.read)) {
  curr.t <- tissues.to.read[i]
  # curr.path.to.read <- paste0(general.path, curr.t, "/count.csv")
  curr.path.to.read <- paste0(general.path, curr.t, "_rm.batch_dge.txt")
  # curr.count <- read.csv(curr.path.to.read, header = T, row.names = 1, stringsAsFactors = F)
  curr.count <- read.table(curr.path.to.read, header = T, row.names = 1, stringsAsFactors = F)
  if (is.null(mca.counts.all.involved)) {
    mca.counts.all.involved <- curr.count
  } else {
    # mca.counts.all.involved <- cbind(mca.counts.all.involved, curr.count)
    mca.counts.all.involved <- merge(mca.counts.all.involved, curr.count, by = "row.names", all = T) 
    rownames(mca.counts.all.involved) <- mca.counts.all.involved$Row.names   
    mca.counts.all.involved <- mca.counts.all.involved[,-1]
    mca.counts.all.involved[is.na(mca.counts.all.involved)] <- 0
  }
}

save(mca.counts.all.involved, file = "mca.counts.all.involved.Rda")
# for (curr.f in curr.path.to.read) {
#   curr.count <- read.table(curr.f, header = T, row.names = 1, stringsAsFactors = F)
#   if (is.null(curr.tissue.count.mtx)) {
#     curr.tissue.count.mtx <- curr.count
#   } else {
#     gene.intersect <- intersect(rownames(curr.tissue.count.mtx), rownames(curr.count))
#     curr.tissue.count.mtx <- cbind(curr.tissue.count.mtx[gene.intersect,], curr.count[gene.intersect,])
#   }
# }



## meta data cleaning
cardiac.rp.all.meta$cell.type.1 <- gsub("\\([^)]*\\)", "", cardiac.rp.all.meta$cell.type)
cardiac.rp.all.meta$cell.type.alone <- unlist(lapply(strsplit(cardiac.rp.all.meta$cell.type.1, "_"), function(x) x[1]))

cardiac.rp.all.meta$cell.type.1 <- tolower(cardiac.rp.all.meta$cell.type.1)
coldata.df <- cardiac.rp.all.meta
save(coldata.df, file = "coldata.df.Rda")

# Construction of a high-resolution reference
ref.list <- construct.high.res.reference(mca.counts.all.involved, coldata.df = coldata.df, criteria = "cell.type.1")
save(ref.list, file = "ref.list.Rda")
# Get expression matrix and meta data of cells used to build the reference, as well as the constructed pseudo-bulk reference
ref.df <- ref.construction(ref.list[[1]], ref.list[[2]], "cell.type")
save(ref.df, file = "ref.df.Rda")

single.round.QP.analysis(ref.df, ref.list[[1]], n.cores = 4, save.to.path = "./", save.to.filename = "stone_et_al_reference_MCA", unix.par = TRUE)
# single.round.QP.analysis(ref.df, stone.et.al, n.cores = 4, save.to.path = "./", save.to.filename = "stone_et_al_test_MCA", unix.par = TRUE)
single.round.QP.analysis(ref.df, stone.et.al.final, n.cores = 4, save.to.path = "./", save.to.filename = "stone_et_al_test_MCA", unix.par = TRUE)


# Read in background and testing identity scores
background.mtx <- read.csv("./stone_et_al_reference_MCA_scale.csv", header = T, row.names = 1, stringsAsFactors = F)
mtx.test <- read.csv("./stone_et_al_test_MCA_scale.csv", header = T, row.names = 1, stringsAsFactors = F)

col.sub <- ncol(background.mtx) - 2

# Conduct reference randomization to get empirical p-value matrix
ref.perc.list <- percentage.calc(background.mtx[,c(1:col.sub)], background.mtx[,c(1:col.sub)])
save(ref.perc.list, file = "ref.perc.list.Rda")

# Conduct test randomization to get empirical p-value matrix
perc.list <- percentage.calc(as.matrix(mtx.test[,c(1:col.sub)]), as.matrix(background.mtx[,c(1:col.sub)]))
save(perc.list, file = "perc.list.Rda")


# Binarization of inference results
bin.count <- binarization.mann.whitney(mtx = mtx.test[,c(1:col.sub)], ref.perc.ls = ref.perc.list, ref.meta = ref.list[[2]], perc.ls = perc.list)
save(bin.count, file = "bin.count.Rda")
# Classificationn
classification <- binary.to.classification(bin.count[,c(1:col.sub)])
rownames(classification) <- classification$barcode
table(classification$call)
save(classification, file = "classification.Rda")
multi.classification.list <- multi.id.curate.qp(binary.counts = bin.count, classification = classification, qp.matrix = mtx.test)
save(multi.classification.list, file = "multi.classification.list.Rda")
# Reassign variables
actual.multi <- multi.classification.list[[1]]
new.classification <- multi.classification.list[[2]]


score.df <- transition.score(actual.multi)
save(score.df, file = "score.df.Rda")


##################################################################################

load("~/cellranger.7.0.1/2.MM/500_5000_10/9.MMsub/3.SCT/6.res0.1_ex_10_12_13/7.MMsub/11.MMsub_copykat/12.GSEA_forloop/c5/MM_celltype.Rda")

p3<-DimPlot(MM,label = T)
p3

score.df$call <- rownames(score.df)
classification |>
  left_join(score.df, by=c("call")) -> classification
rownames(classification) <- classification$barcode

MM <- AddMetaData(MM, classification)
MM$Transition.Score <- MM$entropy
#MM$Transition.Score[is.na(MM$Transition.Score)] <- 0
table(MM$entropy)



features <- c("Transition.Score")
p1 <- FeaturePlot(MM, features = features, cols = c("black", "sienna1"), label = T)
p1
p1 <- FeaturePlot(MM, features = features, cols = c("black", "lightsalmon"), label = T)
p1
p1 <- FeaturePlot(MM, features = features, cols = c("black", "red"), label = F)
p1

pdf("Transition.Score.pdf", width =8, height = 6)
p1
dev.off()


library(ggpubr)
# create a dataset
data <- FetchData(object = MM, vars = c("cellType", "Transition.Score"))

# Plot
data %>%
  ggplot( aes(x=cellType, y=Transition.Score, fill=cellType)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Transition.Score") +
  xlab("") -> p


my_comparisons = list( c("MM0", "MM1"), c("MM0", "MM2"), c("MM0", "MM3"), c("MM0", "MM4"),
                       c("MM1", "MM2"), c("MM1", "MM3"), c("MM1", "MM4"),
                       c("MM2", "MM3"), c("MM2", "MM4"), 
                       c("MM3", "MM4")
                       )
p + stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "t.test") +
    scale_fill_manual(values = c("MM0"="#4b5cc466", 
                                 "MM1"="#16a95166", 
                                 "MM2"="firebrick1", 
                                 "MM3"="firebrick3", 
                                 "MM4"="firebrick4")) -> p2

pdf("Transition.Score.pdf", width = 10, height = 6)
p1 | p2
dev.off()


pdf("Transition.Score3.pdf", width = 20, height = 6)
p3 | p1 | p2
dev.off()

