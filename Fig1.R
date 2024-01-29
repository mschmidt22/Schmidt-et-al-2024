# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig 1

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

## Downloading Data
setwd("")
endocrine_integrated = readRDS('endocrine_integrated.rds')

## Setting Up Object for Integration
DefaultAssay(endocrine_integrated) ="RNA"
split.list <- SplitObject(endocrine_integrated, split.by = "sequencing")
split.list <- lapply(X = split.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = split.list)

## Performing Integration
endocrine.anchors <- FindIntegrationAnchors(object.list = split.list, anchor.features = features)
endocrine_integrated <- IntegrateData(anchorset = endocrine.anchors)
DefaultAssay(endocrine_integrated) <- "integrated"
endocrine_integrated <- ScaleData(endocrine_integrated, vars.to.regress = c('percent.mt', 'nCount_RNA'))
endocrine_integrated <- RunPCA(endocrine_integrated, npcs = 30, verbose = FALSE)
endocrine_integrated <- RunUMAP(endocrine_integrated, reduction = "pca", dims = 1:30)
endocrine_integrated <- FindNeighbors(endocrine_integrated, reduction = "pca", dims = 1:30)
endocrine_integrated <- FindClusters(endocrine_integrated, resolution = 2)

## Assigning Cell Identity
endocrine_integrated = RenameIdents(endocrine_integrated, '0' = "Beta", '1' = "Alpha",
                                    '2' = "Alpha", '3' = "EC", '4' = "Alpha", 
                                    '5' = "Alpha", '6'="EC", '7'="Beta",
                                    '8'="Delta",'9'="Alpha", '10'="Alpha", '11'="Beta",
                                    '12'="Beta",'13' = "PP", '14'="Beta", '15'="Beta",
                                    '16'="Alpha",'17'= "Beta", '18' = "Beta", '19' = "Alpha", '20' = "Beta"
                                    , '21' = "Beta", '22' = "Prolif", '23' = "NeuroEndo", '24' = "Delta"
                                    , '25' = "Alpha", '26' = "Alpha", '27' = "Prog", '28' = "NeuroEndo", 
                                    '29' = "Beta", '30' = "Poly",'31' = "Delta", '32' = "Epsilon", 
                                    '33' = "Alpha", '34' = "Beta", '35' = "Delta", '36' = "Beta", '37' = "NeuroEndo")
endocrine_integrated$endocrine_identity = endocrine_integrated@active.ident

## Identifying Cell-type specific DEGs
DefaultAssay(endocrine_integrated) = "RNA"
Idents(endocrine_integrated) = "endocrine_identity"
DEGs = FindAllMarkers(endocrine_integrated, logfc.threshold = 0.25)

## UMAP (Fig 1B & 1D)
DimPlot(endocrine_integrated, reduction = "umap", repel = TRUE, label = F, cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                                                                                    "#1B9E77", "#7570B3", "#A6761D",
                                                                                    "#E7298A","#FF7F00", "navyblue"))+NoLegend()
DimPlot(endocrine_integrated, reduction = "umap", repel = TRUE, label = F, split.by = "maturation", 
        cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                 "#1B9E77", "#7570B3", "#A6761D", "#E7298A",
                 "#FF7F00", "navyblue"))+NoLegend()

## Splitting by maturation
sc_vitro_cells = subset(endocrine_integrated, maturation %in% "SC-Vitro")
sc_vivo_cells = subset(endocrine_integrated, maturation %in% "SC-Vivo")
adult_cells = subset(endocrine_integrated, maturation %in% "Adult")
fetal_cells = subset(endocrine_integrated, maturation %in% "Fetal")

## Splitting by cell-type (Fig 1C)
sc_vitro_identity = SplitObject(sc_vitro_cells, split.by = "endocrine_identity")
sc_vivo_identity = SplitObject(sc_vivo_cells, split.by = "endocrine_identity")
adult_identity = SplitObject(adult_cells, split.by = "endocrine_identity")
fetal_identity = SplitObject(fetal_cells, split.by = "endocrine_identity")

## Islet Hormone Featureplots (Fig 1E)
DefaultAssay(endocrine_integrated) = "RNA"
FeaturePlot(endocrine_integrated, features = c('INS'))+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(endocrine_integrated, features = c('GCG'))+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(endocrine_integrated, features = c('SST'))+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(endocrine_integrated, features = c('PPY'))+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(endocrine_integrated, features = c('GHRL'))+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(endocrine_integrated, features = c('TPH1'))+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()

## Heatmap of cell-specific markers (Fig 1F)
DefaultAssay(endocrine_integrated) = "RNA"
Idents(endocrine_integrated) <- "endocrine_identity"
endocrine_integrated$endocrine_identity <- factor(endocrine_integrated$endocrine_identity,levels=c('Poly', 'NeuroEndo', 'EC', 'Prolif', 'Prog','Epsilon', 'PP', 'Delta', 'Alpha', 'Beta'))
gene.list <-c('GAP43','RTN1','CNTNAP2','TPH1','FEV','DDC','TOP2A','CENPF','PCNA','NKX2-2','NEUROG3','SOX4','GHRL','ACSL1','PHGR1','PPY','AQP3','ID2', 'SST', 'HHEX', 'LEPR', 'GCG','GC','TTR', 'INS', 'IAPP', 'HADH')
av.exp <- AverageExpression(endocrine_integrated,assays = "RNA",features  = gene.list)
av.exp <- data.matrix(av.exp[[1]], rownames.force = NA)
palette = colorRampPalette(c("gray","white","#DC143C")) (20)
heatmap.2(x = av.exp,dendrogram='none',col= palette, trace="none",density.info="none",key=TRUE,scale="row",Rowv=FALSE,Colv=FALSE,cexCol=1,cexRow = 1)  

## Beta Core Identity Pairwise Comparison (Fig 1G)
DefaultAssay(sc_vitro_cells) = 'RNA'
Idents(sc_vitro_cells) = "endocrine_identity"
sc_vitro_beta_DEGs = FindMarkers(sc_vitro_cells, ident.1 = "Beta", only.pos = T, logfc.threshold = 0.3)
DefaultAssay(sc_vivo_cells) = 'RNA'
Idents(sc_vivo_cells) = "endocrine_identity"
sc_vivo_beta_DEGs = FindMarkers(sc_vivo_cells, ident.1 = "Beta", only.pos = T, logfc.threshold = 0.3)
DefaultAssay(adult_cells) = 'RNA'
Idents(adult_cells) = "endocrine_identity"
adult_beta_DEGs = FindMarkers(adult_cells, ident.1 = "Beta", only.pos = T, logfc.threshold = 0.3)
DefaultAssay(fetal_cells) = 'RNA'
Idents(fetal_cells) = "endocrine_identity"
fetal_beta_DEGs = FindMarkers(fetal_cells, ident.1 = "Beta", only.pos = T, logfc.threshold = 0.3)


