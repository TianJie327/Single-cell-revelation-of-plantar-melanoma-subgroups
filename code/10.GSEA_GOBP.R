library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggraph)
library(magrittr)
library(tibble)
library(tidyverse)
library(SeuratObject)
library(AnnotationHub)
library(BiocGenerics)
library(parallel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(readr)
library(enrichplot)
library(DO.db)
library(dplyr)
library(ggnewscale)
library(ggupset)
library(ggridges)
library(patchwork)
library(tidyverse)
library(biomaRt)
library(DT)

load("MMsub.Rda")

DefaultAssay(object = scRNA) <- "RNA"
melanoma <- scRNA
sample1 <- c("MM0", 
             "MM1", 
             "MM2",
             "MM3",
             "MM4"
      
             )
gsea_GO_BP_symbol_list <- list()

for (i in sample1) {
    markers <- FindMarkers(melanoma, ident.1 = i, min.pct = 0.25, logfc.threshold = 0)
    write.csv(markers, file = paste0(i, "_markers.csv"))
    markers$SYMBOL <- rownames(markers)
    gene <- bitr(rownames(markers), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    GSEA_gene <- gene %>% left_join(markers, by = "SYMBOL")
    rank_list <- GSEA_gene$avg_log2FC
    names(rank_list) <- GSEA_gene$ENTREZID
    rank_list <- sort(rank_list, decreasing = TRUE)
    gsea_GO_BP <- gseGO(rank_list,
        ont           = "BP",
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        nPermSimple   = 10000,
        eps           = 0, 
        minGSSize     = 10,
        maxGSSize     = 500,
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH"
        # seed          = F,
        # verbose       = T,
        # by="fgsea"
    )
    gsea_GO_BP_symbol <- setReadable(gsea_GO_BP, "org.Hs.eg.db", keyType = "ENTREZID")
    gsea_GO_BP_symbol_list[[i]] <- gsea_GO_BP_symbol
    
    save(gsea_GO_BP_symbol, file = paste0(i,"_gsea_GO_BP_symbol.Rda"))
    write.csv(gsea_GO_BP_symbol, file = paste0(i, "_gsea_GO_BP_symbol.csv"), quote = F, row.names = T)
    P_GO_BP <- dotplot(gsea_GO_BP_symbol, showCategory = 20) + ggtitle("dotplot for gsea_GO_BP")
    pdf(file = paste0(i, "_gsea_GO_BP_clusterProfiler.pdf"), width = 10, height = 8)
    print(P_GO_BP)
    dev.off()
    gsea <- as.data.frame(gsea_GO_BP_symbol)
    gsea %>%
        as_tibble() %>%
        arrange(desc(NES)) %>%
        filter(p.adjust < 0.05) %>%
        head(n = 10) -> top5
    top5$group <- i
    group <- top5[, c("Description", "group", "NES", "p.adjust", "setSize")]
    write.csv(group, file = paste0(i, ".csv"), quote = F, row.names = T)
    assign(i, group)
}
save(gsea_GO_BP_symbol_list, file = paste0(i,"_gsea_GO_BP_symbol_list.Rda"))
airqualityL <- rbind(MM0, 
                     MM1, 
                     MM2,
                     MM3,
                     MM4
                     )    
airqualityL$Description <- as.factor(airqualityL$Description) 
airqualityL$Description <- fct_inorder(airqualityL$Description) 
airqualityL <- airqualityL[order(airqualityL$Description, decreasing = T), ] 
airqualityL$Description <- fct_inorder(airqualityL$Description)
save(airqualityL, file = "airqualityL.Rda")
write.csv(airqualityL, file = "airqualityL.csv", quote = F, row.names = T)


load("airqualityL.Rda")
pt.name <- airqualityL
gsea_GO_BP_symbol_list <- list()
load("MM0_gsea_GO_BP_symbol.Rda")
gsea_GO_BP_symbol_list[["MM0"]] <- as.data.frame(gsea_GO_BP_symbol)
load("MM1_gsea_GO_BP_symbol.Rda")
gsea_GO_BP_symbol_list[["MM1"]] <- as.data.frame(gsea_GO_BP_symbol)
load("MM2_gsea_GO_BP_symbol.Rda")
gsea_GO_BP_symbol_list[["MM2"]] <- as.data.frame(gsea_GO_BP_symbol)
load("MM3_gsea_GO_BP_symbol.Rda")
gsea_GO_BP_symbol_list[["MM3"]] <- as.data.frame(gsea_GO_BP_symbol)
load("MM4_gsea_GO_BP_symbol.Rda")
gsea_GO_BP_symbol_list[["MM4"]] <- as.data.frame(gsea_GO_BP_symbol)

MM0 <- gsea_GO_BP_symbol_list[["MM0"]]
MM0$group <- "MM0"
MM1 <- gsea_GO_BP_symbol_list[["MM1"]]
MM1$group <- "MM1"
MM2 <- gsea_GO_BP_symbol_list[["MM2"]]
MM2$group <- "MM2"
MM3 <- gsea_GO_BP_symbol_list[["MM3"]]
MM3$group <- "MM3"
MM4 <- gsea_GO_BP_symbol_list[["MM4"]]
MM4$group <- "MM4"

airqualityL <- rbind(MM0,
                     MM1, 
                     MM2,
                     MM3,
                     MM4
)  


filter <- filter(airqualityL, airqualityL$Description %in% unique(pt.name$Description))
filter$Description <- as.factor(filter$Description)
filter$Description <- fct_inorder(filter$Description) 
filter <- filter[order(filter$Description, decreasing = T), ] 
filter$Description <- fct_inorder(filter$Description)

P3 <- ggplot(filter, aes(x = group, y = Description)) +
    geom_point(aes(
        size = NES,
        color = NES
    )) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "left"
    ) +
    scale_color_gradient(low = "lightblue", high = "red") +
    labs(x = NULL, y = NULL) +
    scale_y_discrete(position = "right")

P3
pdf("GSEA_GOBP.pdf", width = 8)
P3
dev.off()







