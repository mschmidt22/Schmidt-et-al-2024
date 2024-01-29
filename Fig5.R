# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig 5

#Download Packages
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(lubridate)
library(RColorBrewer)
library(enrichR)
library(Matrix)
library(circlize)
library(ComplexHeatmap)
library(gplots)
library(scales)
library(EnhancedVolcano)
library(tidyverse)
library(scclusteval)
library(singleseqgset)
library(heatmap3)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(fmsb)
library(tidyverse)
library(SeuratWrappers)
library(monocle3)
library(data.table)

## Downloading Beta Data
setwd("")
endocrine_integrated = readRDS('endocrine_integrated.rds')
all_beta = subset(endocrine_integrated, endocrine_identity %in% "Beta")
sc_adult_beta = subset(all_beta, maturation %in% c("SC-Vitro", "Adult"))

## Gene set enrichment analysis (Fig 5A)
## https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=GO
h.human <- msigdbr(species="Homo sapiens",category="C5")
h.human = subset(h.human, gs_subcat %in% c('GO:CC'))
h.names <- unique(h.human$gs_name)
h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}
Idents(sc_adult_beta) <- "maturation"
logfc.data <- logFC(cluster.ids=sc_adult_beta@meta.data$maturation,expr.mat=sc_adult_beta@assays[["RNA"]]@data)
names(logfc.data)
gse.res <- wmw_gsea(expr.mat=sc_adult_beta@assays[["RNA"]]@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
names(gse.res)

## Gene set enrichment analysis (Fig 5B)
## https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=GO
h.human <- msigdbr(species="Homo sapiens",category="C5")
h.human = subset(h.human, gs_subcat %in% c('GO:BP'))
h.names <- unique(h.human$gs_name)
h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}
Idents(sc_adult_beta) <- "maturation"
logfc.data <- logFC(cluster.ids=sc_adult_beta@meta.data$maturation,expr.mat=sc_adult_beta@assays[["RNA"]]@data)
names(logfc.data)
gse.res <- wmw_gsea(expr.mat=sc_adult_beta@assays[["RNA"]]@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
names(gse.res)

## Dotplot of Genes involved in Insulin Secretion (Fig 5C)
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
DotPlot(all_beta, features = c('VAMP2', 'VAMP3', 'STX1A', 'STXBP1', 'SNAP25', 'SYT7'), cols = c("lightgrey", "blue"), dot.scale = 10, scale=F)+RotatedAxis()+theme(axis.text.x = element_text(size = 10))

## Neural Identity Violin Panels (Fig 5D)
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
VlnPlot(all_beta, features = c('KIF5C', 'BASP1', 'CD24', 'MARCKS', 'MAP1B', 'MDK',
                               'MLLT11', 'RUFY3', 'NEFM', 'HMGCR', 'NAV1', 'NTM', 'PLPPR5', 'CHODL', 'EVL', 
                               'AKT1', 'ABLIM1', 'GAP43', 'NEFL'), 
        split.by  = "maturation", idents = c("Adult", "SC-Vivo", "SC-Vitro", "Fetal"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = c('H3F3A', 'ISL1', 'GTF2I', 'ONECUT2', 'DPYSL2', 'ASCL2', 'SOX4', 'SOX11', 'RUNX1T1', 'ARID1A',
                               'ARID1B', 'PBX1', 'ELAVL2', 'ASCL1'), 
        split.by  = "maturation", idents = c("Adult", "SC-Vivo", "SC-Vitro", "Fetal"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = c('STMN1', 'STMN2', 'DCX', 'PLXNA2','VLDLR', 'CRMP1', 'LSAMP', 
                               'NCAM1', 'PLXNC1', 'CXCL12','STMN4'), 
        split.by  = "maturation", idents = c("Adult", "SC-Vivo", "SC-Vitro", "Fetal"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = c('MARCKSL1', 'C1QL1', 'RAB3B', 'CADPS', 'CAMK2N1', 'NBEA', 'SYT11',
                               'APOE', 'CALY','CHRNA3'), 
        split.by  = "maturation", idents = c("Adult", "SC-Vivo", "SC-Vitro", "Fetal"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = c('CALB2', 'FXYD6', 'PCP4', 'ATP2B1', 'CELF4', 'CACNA2D1', 'CACNA1A'), 
        split.by  = "maturation", idents = c("Adult", "SC-Vivo", "SC-Vitro", "Fetal"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'))+NoLegend()
