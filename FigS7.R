# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig S7

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
adult_fetal_beta = subset(all_beta, maturation %in% c("Adult", "Fetal"))
sc_cells = subset(endocrine_integrated, maturation %in% c('SC-Vitro', 'SC-Vivo'))
sc_beta_EC = subset(sc_cells, endocrine_identity %in% c("Beta", 'EC'))

## Heatmap of Genes involved in axon, synapse, and dendrite components (Fig S7A)
## Gene Sets from:
## https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOCC_AXON.html
## https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOCC_SYNAPSE.html
## https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOCC_DENDRITIC_TREE.html
axon_genes = read.delim('GOCC_AXON.v2022.1.Hs.txt')
axon_genes = axon_genes$GOCC_AXON[2:644]
synapse_genes = read.delim('GOCC_SYNAPSE.v2022.1.Hs.txt')
synapse_genes = synapse_genes$GOCC_SYNAPSE[2:1343]
dendrite_genes = read.delim('GOCC_DENDRITIC_TREE.v2022.1.Hs.txt')
dendrite_genes = dendrite_genes$GOCC_DENDRITIC_TREE[2:624]
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) <- "maturation"
all_beta$maturation <- factor(all_beta$maturation,levels=c('SC-Vitro', 'SC-Vivo','Adult', 'Fetal'))
av.exp <- AverageExpression(all_beta,assays = "RNA",features  = c(axon_genes, synapse_genes, dendrite_genes))
av.exp <- data.matrix(av.exp[[1]], rownames.force = NA)
palette = colorRampPalette(c("blue","white","red")) (20)
heatmap.2(x = av.exp,dendrogram='none',col= palette, trace="none",density.info="none",key=TRUE,scale="row",Rowv=FALSE,Colv=FALSE,cexCol=1,cexRow = 1)  

## Gene set enrichment analysis (Fig S7B)
## https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=GO
h.human <- msigdbr(species="Homo sapiens",category="C5")
h.human = subset(h.human, gs_subcat %in% c('GO:CC'))
h.names <- unique(h.human$gs_name)
h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}
Idents(adult_fetal_beta) <- "maturation"
logfc.data <- logFC(cluster.ids=adult_fetal_beta@meta.data$maturation,expr.mat=adult_fetal_beta@assays[["RNA"]]@data)
names(logfc.data)
gse.res <- wmw_gsea(expr.mat=adult_fetal_beta@assays[["RNA"]]@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
names(gse.res)

## Gene set enrichment analysis (Fig S7C)
## https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=GO
h.human <- msigdbr(species="Homo sapiens",category="C5")
h.human = subset(h.human, gs_subcat %in% c('GO:BP'))
h.names <- unique(h.human$gs_name)
h.sets <- vector("list",length=length(h.names))
names(h.sets) <- h.names
for (i in names(h.sets)) {
  h.sets[[i]] <- pull(h.human[h.human$gs_name==i,"gene_symbol"])
}
Idents(adult_fetal_beta) <- "maturation"
logfc.data <- logFC(cluster.ids=adult_fetal_beta@meta.data$maturation,expr.mat=adult_fetal_beta@assays[["RNA"]]@data)
names(logfc.data)
gse.res <- wmw_gsea(expr.mat=adult_fetal_beta@assays[["RNA"]]@data,cluster.cells=logfc.data[[1]],log.fc.cluster=logfc.data[[2]],gene.sets=h.sets)
names(gse.res)

## Integration of SC-Beta and SC-EC cells
DefaultAssay(sc_beta_EC) ="RNA"
split.list <- SplitObject(sc_beta_EC, split.by = "sequencing")
split.list <- lapply(X = split.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = split.list)
endocrine.anchors <- FindIntegrationAnchors(object.list = split.list, anchor.features = features)
sc_beta_EC_integrated <- IntegrateData(anchorset = endocrine.anchors)
DefaultAssay(sc_beta_EC_integrated) <- "integrated"
sc_beta_EC_integrated <- ScaleData(sc_beta_EC_integrated, vars.to.regress = c('percent.mt', 'nCount_RNA'))
sc_beta_EC_integrated <- RunPCA(sc_beta_EC_integrated, npcs = 30, verbose = FALSE)
sc_beta_EC_integrated <- RunUMAP(sc_beta_EC_integrated, reduction = "pca", dims = 1:20)
sc_beta_EC_integrated <- FindNeighbors(sc_beta_EC_integrated, reduction = "pca", dims = 1:20)
sc_beta_EC_integrated <- FindClusters(sc_beta_EC_integrated, resolution = 2)

## UMAP (Fig S7D)
DimPlot(sc_beta_EC_integrated, reduction = "umap", repel = TRUE, label = F, group.by = "endocrine_identity", cols = c("#E41A1C", "#4DAF4A"))+NoLegend()

## Violin Plots of Serotonin Associated Genes (Fig S7E)
DefaultAssay(sc_beta_EC_integrated) ="RNA"
Idents(sc_beta_EC_integrated) ="endocrine_identity"
VlnPlot(sc_beta_EC_integrated, features = c('TPH1', 'SLC18A1', 'LMX1A'), pt.size = 0, cols = c("#E41A1C", "#4DAF4A"))

## Downloading Beta & EC Data
all_beta = subset(endocrine_integrated, endocrine_identity %in% "Beta")
all_EC = subset(endocrine_integrated, endocrine_identity %in% "EC")
sc_beta = subset(all_beta, maturation %in% "SC-Vitro")
adult_beta = subset(all_beta, maturation %in% "Adult")
sc_EC = subset(all_EC, maturation %in% "SC-Vitro")
beta_EC = merge(adult_beta, c(sc_beta, sc_EC))

## Addition of sorting Column  into Beta-EC Seurat Object
added_column = c()
length = ncol(beta_EC)
count = 1
while (count <= length) {
  if(beta_EC$endocrine_identity[count] %in% c("EC")){
    value = "SC-EC"
  }
  else if(beta_EC$maturation[count] %in% c("SC-Vitro")){
    value = "SC-Beta"
  }
  else{
    value = "Adult Beta"
  }
  added_column = append(added_column, value)
  count = count + 1
}
View(added_column)
beta_EC$sorting = added_column

## Scaling genes
all.genes <- rownames(beta_EC)
beta_EC <- ScaleData(beta_EC, features = all.genes)

## Neural Identity Violin Panels (Fig S7F)
DefaultAssay(beta_EC) = "RNA"
Idents(beta_EC) = "sorting"
VlnPlot(beta_EC, features = c('KIF5C', 'BASP1', 'CD24', 'MARCKS', 'MAP1B', 'MDK',
                              'MLLT11', 'RUFY3', 'NEFM', 'HMGCR', 'NAV1', 'NTM', 'PLPPR5', 'CHODL', 'EVL', 
                              'AKT1', 'ABLIM1', 'GAP43', 'NEFL'), 
        split.by  = "sorting", idents = c("SC-EC", "SC-Beta", "Adult Beta"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'magenta2', '#4DAF4A'))+NoLegend()
VlnPlot(beta_EC, features = c('H3F3A', 'ISL1', 'GTF2I', 'ONECUT2', 'DPYSL2', 'ASCL2', 'SOX4', 'SOX11', 'RUNX1T1', 'ARID1A',
                              'ARID1B', 'PBX1', 'ELAVL2', 'ASCL1'), 
        split.by  = "sorting", idents = c("SC-EC", "SC-Beta", "Adult Beta"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'magenta2', '#4DAF4A'))+NoLegend()
VlnPlot(beta_EC, features = c('STMN1', 'STMN2', 'DCX', 'PLXNA2','VLDLR', 'CRMP1', 'LSAMP', 
                              'NCAM1', 'PLXNC1', 'CXCL12','STMN4'), 
        split.by  = "sorting", idents = c("SC-EC", "SC-Beta", "Adult Beta"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'magenta2', '#4DAF4A'))+NoLegend()
VlnPlot(beta_EC, features = c('MARCKSL1', 'C1QL1', 'RAB3B', 'CADPS', 'CAMK2N1', 'NBEA', 'SYT11',
                              'APOE', 'CALY','CHRNA3'), 
        split.by  = "sorting", idents = c("SC-EC", "SC-Beta", "Adult Beta"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'magenta2', '#4DAF4A'))+NoLegend()
VlnPlot(beta_EC, features = c('CALB2', 'FXYD6', 'PCP4', 'ATP2B1', 'CELF4', 'CACNA2D1', 'CACNA1A'), 
        split.by  = "sorting", idents = c("SC-EC", "SC-Beta", "Adult Beta"),stack = T, 
        pt.size = 0, cols = c('dodgerblue', 'magenta2', '#4DAF4A'))+NoLegend()