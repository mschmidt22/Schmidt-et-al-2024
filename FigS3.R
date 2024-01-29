# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig S3

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

## Downloading Alpha Data
setwd("")
endocrine_integrated = readRDS('endocrine_integrated.rds')
endocrine_integrated = subset(endocrine_integrated, maturation %in% c("SC-Vitro", "Adult", "Fetal"))
all_alpha = subset(endocrine_integrated, endocrine_identity %in% "Alpha")

## PCA & Clustering
DefaultAssay(all_alpha) = "RNA"
all_alpha <- NormalizeData(all_alpha)
all_alpha <- FindVariableFeatures(all_alpha, selection.method = "vst", nfeatures = 2000)
all_alpha <- ScaleData(all_alpha, vars.to.regress = c('percent.mt', 'nCount_RNA', 'sequencing', 'orig.ident'))
all_alpha <- RunPCA(all_alpha, features = VariableFeatures(object = all_alpha))
all_alpha <- FindNeighbors(all_alpha, dims = 1:20)
all_alpha <- FindClusters(all_alpha, resolution = 1)
all_alpha <- RunUMAP(all_alpha, dims = 1:20)

## UMAP (Fig S3A) 
DimPlot(all_alpha, reduction = "umap", label = F, group.by = "sequencing2",
        cols = c('dodgerblue', '#D95F02', '#7570B3', 'green3', '#E7298A', '#1B9E77'))+NoLegend()

## Pearson Correlation Heatmap (Fig S3B)
Idents(all_alpha) <-"sequencing2"
DefaultAssay(all_alpha) <- 'RNA'
all_alpha$sequencing2 <- factor(all_alpha$sequencing2,levels=c('Fetal', 'Adult', 'Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC'))
all_alpha<-FindVariableFeatures(all_alpha,assay="RNA",nfeatures = 1000) 
genelist <- all_alpha@assays[["RNA"]]@var.features #for most variable genes
av.exp <- AverageExpression(all_alpha,features = genelist,assay="RNA",group.by ="sequencing2")
av.exp <- as.matrix(av.exp[[1]])
cor.exp <- cor(av.exp,method = "pearson")
palette = colorRampPalette(c("steelblue","#f4f0da","darkred")) (200)
heatmap.2(x = cor.exp, col = palette,trace="none",density.info='none',dendrogram="none",Colv=FALSE,Rowv=FALSE) #for printing color key

## Alpha DEGs Heatmap (Fig S3C)
DefaultAssay(all_alpha) = 'RNA'
all_genes = rownames(all_alpha)
all_alpha = ScaleData(all_alpha, features = all_genes)
all_alpha_shuffle = RandomSubsetData(all_alpha, .999)
Idents(all_alpha_shuffle) = "sequencing2"
my_levels <- c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC', 'Fetal', 'Adult')
all_alpha_shuffle@active.ident <- factor(x = all_alpha_shuffle@active.ident, levels = my_levels)
alpha_protocol_DEGs <- FindAllMarkers(all_alpha_shuffle, logfc.threshold = 0.25, only.pos = T)
alpha_protocol_DEGs %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
my_levels <- c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC', 'Fetal', 'Adult')
all_alpha_shuffle@active.ident <- factor(x = all_alpha_shuffle@active.ident, levels = my_levels)
DoHeatmap(all_alpha_shuffle, features = top20$gene, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A"))

## Violin Plots of Alpha markers (Fig S3D)
DefaultAssay(all_alpha) = "RNA"
Idents(all_alpha) = "sequencing2"
my_levels <- c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC', 'Fetal', 'Adult')
all_alpha@active.ident <- factor(x = all_alpha@active.ident, levels = my_levels)
VlnPlot(all_alpha, features = c("GCG"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#1B9E77', 'green3', 'dodgerblue'),pt.size = 0)+NoLegend()
VlnPlot(all_alpha, features = c("TTR"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#1B9E77', 'green3', 'dodgerblue'),pt.size = 0)+NoLegend()
VlnPlot(all_alpha, features = c("IRX2"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#1B9E77', 'green3', 'dodgerblue'),pt.size = 0)+NoLegend()
VlnPlot(all_alpha, features = c("ARX"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#1B9E77', 'green3', 'dodgerblue'),pt.size = 0)+NoLegend()

## Downloading Delta Data
setwd("")
endocrine_integrated = readRDS('endocrine_integrated.rds')
endocrine_integrated = subset(endocrine_integrated, maturation %in% c("SC-Vitro", "Adult", "Fetal"))
all_delta = subset(endocrine_integrated, endocrine_identity %in% "Delta")

## PCA & Clustering
DefaultAssay(all_delta) = "RNA"
all_delta <- NormalizeData(all_delta)
all_delta <- FindVariableFeatures(all_delta, selection.method = "vst", nfeatures = 2000)
all_delta <- ScaleData(all_delta, vars.to.regress = c('percent.mt', 'nCount_RNA', 'sequencing', 'orig.ident'))
all_delta <- RunPCA(all_delta, features = VariableFeatures(object = all_delta))
all_delta <- FindNeighbors(all_delta, dims = 1:20)
all_delta <- FindClusters(all_delta, resolution = 1)
all_delta <- RunUMAP(all_delta, dims = 1:20)

## UMAP (Fig S3E) 
DimPlot(all_delta, reduction = "umap", label = F, group.by = "sequencing2",
        cols = c('dodgerblue', '#D95F02', '#7570B3', 'green3', '#E7298A', '#1B9E77'))+NoLegend()

## Pearson Correlation (Fig S3F)
Idents(all_delta) <-"sequencing2"
DefaultAssay(all_delta) <- 'RNA'
all_delta$sequencing2 <- factor(all_delta$sequencing2,levels=c('Fetal', 'Adult', 'Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC'))
all_delta<-FindVariableFeatures(all_delta,assay="RNA",nfeatures = 1000) 
genelist <- all_delta@assays[["RNA"]]@var.features #for most variable genes
av.exp <- AverageExpression(all_delta,features = genelist,assay="RNA",group.by ="sequencing2")
av.exp <- as.matrix(av.exp[[1]])
cor.exp <- cor(av.exp,method = "pearson")
palette = colorRampPalette(c("steelblue","#f4f0da","darkred")) (200)
heatmap.2(x = cor.exp, col = palette,trace="none",density.info='none',dendrogram="none",Colv=FALSE,Rowv=FALSE) #for printing color key

## Delta DEGs Heatmap (Fig S3G)
DefaultAssay(all_delta) = 'RNA'
all_genes = rownames(all_delta)
all_delta = ScaleData(all_delta, features = all_genes)
all_delta_shuffle = RandomSubsetData(all_delta, .999)
Idents(all_delta_shuffle) = "sequencing2"
my_levels <- c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC', 'Fetal', 'Adult')
all_delta_shuffle@active.ident <- factor(x = all_delta_shuffle@active.ident, levels = my_levels)
delta_protocol_DEGs <- FindAllMarkers(all_delta_shuffle, logfc.threshold = 0.25, only.pos = T)
delta_protocol_DEGs %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
my_levels <- c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC', 'Fetal', 'Adult')
all_delta_shuffle@active.ident <- factor(x = all_delta_shuffle@active.ident, levels = my_levels)
DoHeatmap(all_delta_shuffle, features = top20$gene, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A"))



## Violin Plots of Delta markers (Fig S3H)
DefaultAssay(all_delta) = "RNA"
Idents(all_delta) = "sequencing2"
my_levels <- c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC', 'Fetal', 'Adult')
all_delta@active.ident <- factor(x = all_delta@active.ident, levels = my_levels)
VlnPlot(all_delta, features = c("SST"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#1B9E77', 'green3', 'dodgerblue'),pt.size = 0)+NoLegend()
VlnPlot(all_delta, features = c("RBP4"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#1B9E77', 'green3', 'dodgerblue'),pt.size = 0)+NoLegend()
VlnPlot(all_delta, features = c("HHEX"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#1B9E77', 'green3', 'dodgerblue'),pt.size = 0)+NoLegend()
VlnPlot(all_delta, features = c("SEC11C"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#1B9E77', 'green3', 'dodgerblue'),pt.size = 0)+NoLegend()


