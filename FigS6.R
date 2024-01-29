# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig S6

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
sc_vitro_vivo_beta = subset(all_beta, maturation %in% c('SC-Vitro', 'SC-Vivo'))

## Lineplots of Expression level of TFs before and after transplantation (Fig S6A)
DefaultAssay(sc_vitro_vivo_beta) = 'RNA'
all_genes = rownames(sc_vitro_vivo_beta)
sc_vitro_vivo_beta = ScaleData(sc_vitro_vivo_beta, features = all_genes)
DefaultAssay(sc_vitro_vivo_beta) = "RNA"
Idents(sc_vitro_vivo_beta) = "sequencing"
DE_TFs = as.data.frame(AverageExpression(sc_vitro_vivo_beta, features = c('CTNNB1', 'FOXA1', 'FOXA2',
                                                                     'ISL1', 'NKX2-2', 'NKX6-1',
                                                                     'ONECUT2', 'PAX4', 'PBX1',
                                                                     'PDX1', 'PROX1', 'SOX4'), assays = "RNA", slot = "scale.data"))

## Regulon Analysis (Fig SB-D)

### Setting up regulon 
build_loom(file.name = "input.loom",
           dgem = sc_vitro_vivo_beta@assays$RNA@counts,
           title = "SC-Vitro & SC-Vivo Beta Comparison",
           default.embedding = sc_vitro_vivo_beta@reductions$umap@cell.embeddings,
           default.embedding.name = "UMAP")
loom <- open_loom("input.loom", mode = "r+")
add_col_attr(loom = loom, key = "orig.ident", value = sc_vitro_vivo_beta@meta.data$orig.ident)
add_col_attr(loom = loom, key = "cell_label", value = sc_vitro_vivo_beta@meta.data$cell_label)
add_col_attr(loom = loom, key = "age", value = sc_vitro_vivo_beta@meta.data$age)
add_col_attr(loom = loom, key = "sequencing", value = sc_vitro_vivo_beta@meta.data$sequencing)
add_col_attr(loom = loom, key = "maturation", value = sc_vitro_vivo_beta@meta.data$maturation)
add_col_attr(loom = loom, key = "endocrine_identity", value = sc_vitro_vivo_beta@meta.data$endocrine_identity)
close_loom(loom)

### SCENIC regulon run in command line interface
### According to Aibar et al. 2017 (DOI: 10.1038/nmeth.4463)

### Read Info from .loom file
output.loom <- open_loom('output.loom')
exprMat <- get_dgem(output.loom)
regulons_incidMat <- get_regulons(output.loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(output.loom, column.attr.name='RegulonsAUC')
embeddings <- get_embeddings(output.loom)
orig.ident <- get_cell_annotation(output.loom,annotations.columns = "orig.ident")
cell_label <- get_cell_annotation(output.loom,annotations.columns = "cell_label")
age <- get_cell_annotation(output.loom,annotations.columns = "age")
sequencing <- get_cell_annotation(output.loom,annotations.columns = "sequencing")
maturation <- get_cell_annotation(output.loom,annotations.columns = "maturation")
endocrine_identity <- get_cell_annotation(output.loom,annotations.columns = "endocrine_identity")
close_loom(output.loom)

### Add Regulon Data to Seurat Object
AUCmat <- AUCell::getAUC(regulonAUC)
sc_vitro_vivo_beta[["AUC"]] <- CreateAssayObject(data = AUCmat)
DefaultAssay(sc_vitro_vivo_beta) <- "AUC"
sc_vitro_vivo_beta <- ScaleData(sc_vitro_vivo_beta, assay = "AUC")

### Calculate & Plotting RSS S-plots (Fig S6B)
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=maturation[colnames(regulonAUC), "maturation"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
options(ggrepel.max.overlaps = Inf)
plotRSS_oneSet <- function(rss, setName, n=15)
{
  library(ggplot2)
  library(ggrepel)
  
  rssThisType <- sort(rss[,setName], decreasing=TRUE)
  thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
  thisRss$regulon[(n+1):nrow(thisRss)] <- NA
  
  ggplot(thisRss, aes(x=rank, y=rss)) + 
    ggtitle(setName) + 
    geom_point(color='magenta2',size=3)+
    geom_label_repel(size = 5,aes(label = regulon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'black',
                     na.rm=TRUE,
                     segment.alpha=1,
                     alpha=1,
                     seed = 555)+
    theme_classic()
}
plotRSS_oneSet(rss, setName = "SC-Vitro", n =10)
plotRSS_oneSet <- function(rss, setName, n=15)
{
  library(ggplot2)
  library(ggrepel)
  
  rssThisType <- sort(rss[,setName], decreasing=TRUE)
  thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
  thisRss$regulon[(n+1):nrow(thisRss)] <- NA
  
  ggplot(thisRss, aes(x=rank, y=rss)) + 
    ggtitle(setName) + 
    geom_point(color='darkorchid4',size=3)+
    geom_label_repel(size = 5,aes(label = regulon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'black',
                     na.rm=TRUE,
                     segment.alpha=1,
                     alpha=1,
                     seed = 555)+
    theme_classic()
}
plotRSS_oneSet(rss, setName = "SC-Vivo", n =10)


### Regulon Activity Score before and after Transplantation (Fig S6C)
selectedResolution <- 'maturation'
cells_per_maturation <- split(rownames(maturation), maturation[,selectedResolution]) 
regulonActivity_by_maturation <- sapply(cells_per_maturation,
                                        function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

### Regulon UMAP (Fig S6D)
DefaultAssay(sc_vitro_vivo_beta) ="RNA"
sc_vitro_vivo_beta <- NormalizeData(sc_vitro_vivo_beta, normalization.method = "LogNormalize", scale.factor = 10000)
sc_vitro_vivo_beta <- FindVariableFeatures(sc_vitro_vivo_beta, selection.method = "vst", nfeatures = 2000)
sc_vitro_vivo_beta <- ScaleData(sc_vitro_vivo_beta, vars.to.regress = c('percent.mt', 'nCount_RNA', 'orig.ident'))
sc_vitro_vivo_beta <- RunPCA(sc_vitro_vivo_beta, features = VariableFeatures(object = sc_vitro_vivo_beta))
sc_vitro_vivo_beta <- FindNeighbors(sc_vitro_vivo_beta, dims = 1:18)
sc_vitro_vivo_beta <- FindClusters(sc_vitro_vivo_beta, resolution = 1)
sc_vitro_vivo_beta <- RunUMAP(sc_vitro_vivo_beta, dims = 1:18)
DimPlot(sc_vitro_vivo_beta, reduction = "umap", label = F, group.by = 'maturation', cols = c('magenta2', 'darkorchid4'))+NoLegend()

### Expression & Activity Featureplots (Fig S6D) 
DefaultAssay(sc_vitro_vivo_beta) = "RNA"
FeaturePlot(sc_vitro_vivo_beta, features = "SOX4")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(sc_vitro_vivo_beta, features = "SOX11")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(sc_vitro_vivo_beta, features = "PBX1")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(sc_vitro_vivo_beta, features = "PAX4")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
DefaultAssay(sc_vitro_vivo_beta) = "AUC"
FeaturePlot(sc_vitro_vivo_beta, features = "SOX4(+)")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+NoLegend()
FeaturePlot(sc_vitro_vivo_beta, features = "SOX11(+)")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+NoLegend()
FeaturePlot(sc_vitro_vivo_beta, features = "PBX1(+)")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+NoLegend()
FeaturePlot(sc_vitro_vivo_beta, features = "PAX4(+)")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+NoLegend()

