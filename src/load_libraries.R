suppressMessages(library(DESeq2))
suppressMessages(library(pheatmap))
suppressMessages(library(gplots))
suppressMessages(library(UpSetR))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(plyr))
require(graphics)
require(grDevices)
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(rentrez))
suppressMessages(library(goseq))
suppressMessages(library(dendextend))
suppressMessages(library(GOSemSim))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(GO.db))
suppressMessages(library(annotate))
suppressMessages(library(moe430a.db))
suppressMessages(library(WGCNA))
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
suppressMessages(library(igraph))
#suppressMessages(library(clusterProfiler))
suppressMessages(library(threejs))
suppressMessages(library(htmlwidgets))
suppressMessages(library(pathview))
suppressMessages(library(xlsx))
suppressMessages(library(limma))
suppressMessages(library(reshape))
suppressMessages(library(plotly))
#suppressMessages(library(RamiGO))
suppressMessages(library(ggfortify))
suppressMessages(library(gridExtra))
suppressMessages(library(sinaplot))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(factoextra))
suppressMessages(library(tidyr))
suppressMessages(library(ggrepel))
suppressMessages(library(GSVA))
suppressMessages(library(scales))
suppressMessages(library(dorothea))
MF_GO = godata('org.Mm.eg.db', ont="MF")
CC_GO = godata('org.Mm.eg.db', ont="CC")
BP_GO = godata('org.Mm.eg.db', ont="BP")