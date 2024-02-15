setwd("~/cellranger.7.0.1/2.MM/500_5000_10/20.article_figs/11.supplementary/data")
.libPaths(c("/gdata01/apps/miniconda3/envs/R.4.1/lib/R/library","/gdata01/user/tianjie/R/x86_64-conda-linux-gnu-library/4.1"))

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
library(stringr)
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
library(gridExtra)
library(pheatmap)

# fig1.b ----

load("MM.Rda")

color = c(
  "#db3a34",  # "Endothelial",            #0
  "#131200",  # "Melanoma",                        #1
  "#009ffd",  # "T Cells" ,                  #2
  "#893168",  # "Pericytes",                     #3
  "#9d9171",  # "Keratinocytes",               #4
  "#00a676", # "Fibroblasts",            #5
  "#8a4fff", # "Monocytes and acrophages",                #6
  "#388659", # "Mast Cells",             #7
  "#593f62", # "Neutrophils", #8
  "#f4743b", # "Epidermal appendages keratinocytes", #9
  "#f0d3f7",  # "Lymphatic Endothelial Cells", #10
  "#bd93d8" # "B Cells" #11
) 

fig1b <- DimPlot(MM, reduction = "umap", 
              label = FALSE, 
              pt.size = 0.5 , 
              label.size = 6, 
              raster=FALSE,
              cols = color
) 
fig1b
pdf(file = "fig1b.pdf", width = 11)
fig1b
dev.off()
save(MM, file = "MM.Rda")

# fig1.c ----
load("~/cellranger.7.0.1/2.MM/500_5000_10/20.article_figs/4.fig1/NP_all.celltype.article.Rda")
NP <- NP_all

color <- c(  
  "#db3a34",  #          "Endothelial",            #0
  "#00a676",  #          "Fibroblasts",                        #1
  "#9d9171",  #          "Keratinocytes" ,                  #2
  "#893168",  #          "Pericytes",                     #3
  "#f4743b",  #          "Epidermal appendages keratinocytes",               #4
  "#f0d3f7",  #          "Lymphatic Endothelial Cells",            #5
  "#8a4fff",  #          "Monocytes and Macrophages",                #6
  "#131200",  #          "Melanocytes",             #7
  "#009ffd",  #          "T Cells", #8
  "#388659",  #          "Mast Cells", #9
  "#6BD425"   #          "Schwann Cells" #10
)

fig1c <- DimPlot(NP_all, reduction = "umap", 
                 pt.size = 0.5 , 
                 label.size = 6, 
                 repel = T, 
                 cols = color)

fig1c
pdf(file = "fig1c.pdf", width = 11)
fig1c
dev.off()
save(NP, file = "NP.Rda")

# fig1.d ----

load("fig1d.data.mm.np.Rda")
save(data, file = "fig1d.data.mm.np.Rda")
write.csv(data, file = "fig1d.data.mm.np.csv")
colnames(data)
data <- data[order(data$group, decreasing = F), ]
data$sample <- factor(data$sample, levels = data$sample)

data |>
  pivot_longer(`Endothelial`:`Schwann Cells`, 
               names_to = "celltype", 
               values_to = "count") -> data

data$count <- as.numeric(data$count)

celltype.levels <- c(
  "Endothelial",            #0
  "Fibroblasts",                        #1
  "Keratinocytes" ,                  #2
  "Pericytes",                     #3
  "Epidermal Appendages Keratinocytes",               #4
  "Lymphatic Endothelial Cells",            #5
  "Monocytes and Macrophages",                #6
  "Melanoma or Melanocytes",             #7
  "T Cells", #8
  "Mast Cells", #9
  "Schwann Cells", #10
  "Neutrophils",
  "B Cells"
)

color3 = c(
  "#db3a34",  #          "Endothelial",            #0
  "#00a676",  #          "Fibroblasts",                        #1
  "#9d9171",  #          "Keratinocytes" ,                  #2
  "#893168",  #          "Pericytes",                     #3
  "#f4743b",  #          "Epidermal appendages keratinocytes",               #4
  "#f0d3f7",  #          "Lymphatic Endothelial Cells",            #5
  "#8a4fff",  #          "Monocytes and Macrophages",                #6
  "#131200",  #          "Melanocytes",             #7
  "#009ffd",  #          "T Cells", #8
  "#388659",  #          "Mast Cells", #9
  "#6BD425",   #          "Schwann Cells" #10
  "#593f62", #            "Neutrophils", #8
  "#bd93d8" #             "B Cells" #11
) 

celltype.color <- data.frame(celltype.levels, color3)
row.names(celltype.color) <- celltype.color$celltype.levels

celltype.color <- celltype.color[c(
  "Schwann Cells", 
  "T Cells",  
  "Monocytes and Macrophages", 
  "Mast Cells",
  "B Cells",
  "Neutrophils",
  "Epidermal Appendages Keratinocytes", 
  "Keratinocytes",               
  "Pericytes",                     
  "Fibroblasts",     
  "Lymphatic Endothelial Cells", 
  "Endothelial",   
  "Melanoma or Melanocytes"                       
),]
celltype.color
data$celltype <- factor(data$celltype, levels=celltype.color$celltype.levels)
levels(data$celltype)
mycolors <- celltype.color$color3

fig1.d1 <- ggplot(data=data,aes(x=sample,y=count,fill=celltype)) + 
  geom_bar(stat="identity",position="fill") + 
  scale_fill_manual(values=mycolors)+
  scale_y_continuous(expand = expansion(mult=c(0.01,0.02)),
                     labels = scales::percent_format())+
  labs(
    y="Relative Abundance",
    fill=" ",title="")+
  theme_minimal()+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=10))+
  theme(axis.title.x=element_text(size=14))+
  theme(axis.text.x = element_text(size=12, color = "black"))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.ticks.length=unit(0.3,"cm"))+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11))+
  theme(legend.position = "none")

fig1.d2 <- ggplot(data=data,aes(x=sample,y=count,fill=celltype)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=mycolors)+
  labs(x="Samples",y="Cell Count",fill=" ",title="") +
  theme_minimal()+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(legend.text=element_text(size=10))+
  theme(axis.ticks.length=unit(0.3,"cm"))+
  theme(legend.position = "top",
        legend.spacing.x=unit(1, 'cm'),
        legend.spacing.y=unit(1, 'cm')
  )

fig1d <- fig1.d2/fig1.d1
fig1d
pdf("fig1d.pdf", width = 15, height = 10 )
fig1d
dev.off()

# fig1.e ----
load("fig1.e.group3.vs.Rda")
save(data, file = "fig1.e.group3.vs.Rda")
write.csv(data, file = "fig1.e.group3.vs.csv")
celltype.levels <- c(
  "Endothelial",            #0
  "Fibroblasts",                        #1
  "Keratinocytes" ,                  #2
  "Pericytes",                     #3
  "Epidermal Appendages Keratinocytes",               #4
  "Lymphatic Endothelial Cells",            #5
  "Monocytes and Macrophages",                #6
  "Melanoma or Melanocytes",             #7
  "T Cells", #8
  "Mast Cells", #9
  "Schwann Cells", #10
  "Neutrophils",
  "B Cells"
)


color3 = c(
  "#db3a34",  #          "Endothelial",            #0
  "#00a676",  #          "Fibroblasts",                        #1
  "#9d9171",  #          "Keratinocytes" ,                  #2
  "#893168",  #          "Pericytes",                     #3
  "#f4743b",  #          "Epidermal appendages keratinocytes",               #4
  "#f0d3f7",  #          "Lymphatic Endothelial Cells",            #5
  "#8a4fff",  #          "Monocytes and Macrophages",                #6
  "#131200",  #          "Melanocytes",             #7
  "#009ffd",  #          "T Cells", #8
  "#388659",  #          "Mast Cells", #9
  "#6BD425",   #          "Schwann Cells" #10
  "#593f62", # "Neutrophils", #8
  "#bd93d8" # "B Cells" #11
  
) 


celltype.color <- data.frame(celltype.levels, color3)
row.names(celltype.color) <- celltype.color$celltype.levels
str(celltype.color)

celltype <- celltype.color$celltype.levels
fill.colors <- celltype.color$color3
str(data)

categories <- c("In Situ", "Invasive", "Normal")
my_comparisons <- list()
for (i in 1:(length(categories) - 1)) {
  for (j in (i + 1):length(categories)) {
    my_comparisons <- c(my_comparisons, list(c(categories[i], categories[j])))
  }
}


celltype.box <- list()
data$group3 <- factor(data$group3, levels = c("Normal", "In Situ", "Invasive"))
for (i in 1:length(celltype)) { 
  
  # Create Data
  cell.name <-  celltype[i]
  data.box <- filter(data, celltype == cell.name)
  
  #Basic piechart
  p <- ggboxplot(data.box, x = "group3", y = "percent", color = "black", add = "jitter", fill = fill.colors[i]) +
    xlab("") +
    ylab("Proportion of cells (%)") +
    ggtitle(celltype[i]) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none") + 
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "t.test")
  
  
  celltype.box[[i]] <- print(p)
}

data.box <- filter(data, celltype == "B Cells")
p <- ggboxplot(data.box, x = "group3", y = "percent", color = "black", add = "jitter", fill = "#bd93d8") +
  xlab("") +
  ylab("Proportion of cells (%)") +
  ggtitle("B Cells") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") + 
  stat_compare_means(comparisons = list(c("In Situ", "Invasive")),
                     label = "p.signif",
                     method = "t.test")
p
celltype.box[[13]] <- p

data.box <- filter(data, celltype == "Neutrophils")
p <- ggboxplot(data.box, x = "group3", y = "percent", color = "black", add = "jitter", fill = "#593f62") +
  xlab("") +
  ylab("Proportion of cells (%)") +
  ggtitle("Neutrophils") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") + 
  stat_compare_means(comparisons = list(c("In Situ", "Invasive")),
                     label = "p.signif",
                     method = "t.test")
p
celltype.box[[12]] <- p

P1 <- ggarrange(celltype.box[[8]],celltype.box[[7]],celltype.box[[9]],
                celltype.box[[4]],celltype.box[[1]],celltype.box[[10]],
                celltype.box[[2]],celltype.box[[3]],celltype.box[[5]],
                ncol=3, nrow = 3)
P2 <- ggarrange(celltype.box[[6]],celltype.box[[13]],celltype.box[[12]],celltype.box[[11]], ncol = 4)
fig1.e <- ggarrange(P1, P2, ncol = 1, heights = c(3,1))
fig1.e

pdf("fig1.e.pdf", width = 20, height = 15)
fig1.e
dev.off()


# fig2.b ----
load("fig2.b.webcsea.Rda")
pheatmap(data,
         scale = "row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(255), 
         border_color = NA, 
         filename = "heatmap.pdf",
         width = 8,
         height = 5
)

# fig2.c ----

load("fig2.c.data.Rda")
data <- data[order(data$Breslow, decreasing = F), ]
data$sample <- factor(data$sample, levels = data$sample)
head(data)

data |>
  pivot_longer(`Diploid`:`MM4`, 
               names_to = "celltype", 
               values_to = "count") -> data

data$count <- as.numeric(data$count)
data$celltype <- as.factor(data$celltype)
levels(data$celltype)
data$sample <- as.factor(data$sample)

mycolors <- c(
  "#b4b8ab",
  "#92140C",
  "#000000",
  "#FF8200",
  "#8CD867",
  "#A09BE7")

#PA
PA <- ggplot(data=data,aes(x=sample,y=count,fill=celltype)) + 
  geom_bar(stat="identity",position="fill") + 
  scale_fill_manual(values=mycolors)+
  scale_y_continuous(expand = expansion(mult=c(0.01,0.02)),
                     labels = scales::percent_format())+
  labs(
    y="Relative Abundance",
    fill=" ",title="")+
  theme_minimal()+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(legend.text=element_text(size=10))+
  theme(axis.ticks.length=unit(0.3,"cm"))
PA


PB <- ggplot(data=data,aes(x=sample,y=count,fill=celltype)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values=mycolors)+
  labs(x="",y="Cell Count",fill=" ",title="") +
  theme_minimal()+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=10))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.ticks.length=unit(0.3,"cm"))
PB

load("fig2.c.data.Rda")
data <- data[order(data$Breslow, decreasing = F), ]
data$sample <- factor(data$sample, levels = data$sample)

PC <- ggplot(data=data,aes(x=sample,y=Breslow, fill=Invasive)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("#3bc14a","#8367c7"), labels = c("In Situ","Invasion"))+
  labs(x="Sample",y="Breslow(mm)",fill=" ",title="") +
  theme_minimal()+
  theme(axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=10))+
  theme(axis.text.x = element_text(size = 12, color = "black"))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.ticks.length=unit(0.3,"cm"))+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11))
PC

pdf("fig2.c.pdf", width = 10, height = 10)
PB/PA/PC
dev.off()


# fig4 pie ----

load("fig4.out.Rda")
data <- list()
regions <- colnames(fig4.out)[1:4]

for (i in regions) {
  
  p <- ggplot(fig4.out, aes(x="", y=!!sym(i), fill=Group_outer)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = fig4.out$Color_outer)  +
    theme_void() +
    ggtitle(paste0(i)) + theme(plot.title = element_text(hjust = 0.5))
  print(p)
  
  data[[i]] <- print(p)
}


load("fig4.in.Rda")
data2 <- list()
regions <- colnames(fig4.in)[2:5]
for (i in regions) {
  
  p <- ggplot(fig4.in, aes(x="", y=!!sym(i), fill=Group_inter)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = fig4.in$Colors_inter)  +
    theme_void() +
    ggtitle(paste0(i)) + theme(plot.title = element_text(hjust = 0.5))
  print(p)
  
  data2[[i]] <- print(p)
}


p1 <- plot_grid(plotlist = data)
p2 <- plot_grid(plotlist = data2)

pdf("fig4.outer.pdf", height = 20, width = 20)
p1
dev.off()

pdf("fig4.inter.pdf", height = 10, width = 10)
p2
dev.off()


# fig5.a ----

load("fig5.a.out.Rda")
data <- list()
regions <- colnames(fig5.a.out)[1:4]
for (i in regions) {
  
  p <- ggplot(fig5.a.out, aes(x="", y=!!sym(i), fill=Group_outer)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = fig5.a.out$Color_outer)  +
    theme_void() +
    ggtitle(paste0(i)) + theme(plot.title = element_text(hjust = 0.5))
  print(p)
  
  data[[i]] <- print(p)
}


plot_grid(plotlist = data)

load("fig5.a.in.Rda")
data2 <- list()
regions <- colnames(fig5.a.in)[2:5]
for (i in regions) {
  
  p <- ggplot(fig5.a.in, aes(x="", y=!!sym(i), fill=Group_inter)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = fig5.a.in$Color_inter)  +
    theme_void() +
    ggtitle(paste0(i)) + theme(plot.title = element_text(hjust = 0.5))
  print(p)
  
  data2[[i]] <- print(p)
}


plot_grid(plotlist = data2)

p1 <- plot_grid(plotlist = data)
p2 <- plot_grid(plotlist = data2)

pdf("fig5.a.outer.pdf", height = 20, width = 20)
p1
dev.off()

pdf("fig5.a.inter.pdf", height = 10, width = 10)
p2
dev.off()



# fig5.b ----
load("fig5.b.data.box.Rda")
my_comparisons <- list()
categories <- c("Normal", "In Situ", "Invasive", "Metastatic")
for (i in 1:(length(categories) - 1)) {
  for (j in (i + 1):length(categories)) {
    my_comparisons <- c(my_comparisons, list(c(categories[i], categories[j])))
  }
}
colnames(data.box)

compar.pt <- list()

for (i in colnames(data.box)[21:24]) {
  
  compar.pt[[i]] <- ggboxplot(data.box, x = "Category", y = i, color = "Category", add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,
                       label = "p.signif",
                       method = "t.test") +
    xlab("") +
    ylab("Proportion of cells") +
    ggtitle(i) + # 添加标题
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "none") 
  
}

#library(cowplot)
plot_grid(plotlist = compar.pt)



# fig5.c ----
load("fig5.cd.data.box.Rda")

time_col <- unlist(data.box[,paste0("PFStime")])
event_col <- unlist(data.box[, paste0("PFS")])


km.fit <- survfit(Surv(time_col, event_col) ~ Max_Column_Name, data=data.box, se.fit = TRUE, conf.int=.95, conf.type='log-log')
summary(km.fit)
labs <- names(table(data.box$Max_Column_Name))
p0 <- ggsurvplot(km.fit, 
                 conf.int = TRUE, 
                 palette = "aaas", 
                 risk.table = TRUE, 
                 risk.table.col="strata", 
                 pval=TRUE, 
                 fun=function(y) y*100,
                 legend.labs = labs,
                 xlab = "Time in days",
                 ylab = "Progress Free Survival(%)",
                 ggtheme = theme_bw()
)
p0

pdf(paste0("fig5.d_PFS.pdf"), height  = 6, width = 6)
print(p0)
dev.off()


anti <- "MITF_CXCL1_ECRG4_NGFR"
cutoff = data.box[,anti] <= 0.0103092783505155
time_col <- unlist(data.box[,paste0("OStime")])
event_col <- unlist(data.box[, paste0("OS")])

fit <- survfit(Surv(time_col, event_col) ~cutoff, data = td)

p1 <- ggsurvplot(
  fit, 
  conf.int = TRUE, 
  palette = "aaas", 
  risk.table = TRUE, 
  risk.table.col="strata", 
  pval=TRUE, 
  fun=function(y) y*100,
  legend.labs = c("High", "Low"),
  xlab = "Time in days",
  ylab = "Overall Survival(%)",
  ggtheme = theme_bw()
)
p1
pdf(paste0("fig5.c_OS_MITF_CXCL1_ECRG4_NGFR.pdf"), height  = 6, width = 6)
print(p1)
dev.off()



