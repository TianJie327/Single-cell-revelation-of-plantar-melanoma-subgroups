library(SoupX)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggraph)
library(limma)
library(magrittr)
library(tibble)

#一、SoupX
#==========================================
setwd("/gdata01/user/tianjie/cellranger.7.0.1/soupx/MM2/MM2/outs")

pt = "/gdata01/user/tianjie/cellranger.7.0.1/soupx/MM2/MM2/outs"
sc = load10X(pt)


sc = setContaminationFraction(sc, 0.2)
out = adjustCounts(sc)
DropletUtils:::write10xCounts("./MM2_strainedCounts_0.2", out)



