# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig 2

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
endocrine_integrated = subset(endocrine_integrated, maturation %in% c("SC-Vitro", "Adult", "Fetal"))
all_beta = subset(endocrine_integrated, endocrine_identity %in% "Beta")

## PCA & Clustering
DefaultAssay(all_beta) = "RNA"
all_beta <- NormalizeData(all_beta)
all_beta <- FindVariableFeatures(all_beta, selection.method = "vst", nfeatures = 2000)
all_beta <- ScaleData(all_beta, vars.to.regress = c('percent.mt', 'nCount_RNA', 'sequencing', 'orig.ident'))
all_beta <- RunPCA(all_beta, features = VariableFeatures(object = all_beta))
all_beta <- FindNeighbors(all_beta, dims = 1:20)
all_beta <- FindClusters(all_beta, resolution = 1)
all_beta <- RunUMAP(all_beta, dims = 1:20)

## UMAP (Fig 2A) 
DimPlot(all_beta, reduction = "umap", label = F, group.by = "sequencing2",
        cols = c('dodgerblue', '#D95F02', '#7570B3', 'green3', '#E7298A', '#1B9E77'))+NoLegend()

## Pearson Correlation Heatmap (Fig 2B)
Idents(all_beta) <-"sequencing2"
DefaultAssay(all_beta) <- 'RNA'
all_beta$sequencing2 <- factor(all_beta$sequencing2,levels=c('Fetal', 'Adult', 'Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC'))
all_beta<-FindVariableFeatures(all_beta,assay="RNA",nfeatures = 1000) 
genelist <- all_beta@assays[["RNA"]]@var.features #for most variable genes
av.exp <- AverageExpression(all_beta,features = genelist,assay="RNA",group.by ="sequencing2")
av.exp <- as.matrix(av.exp[[1]])
cor.exp <- cor(av.exp,method = "pearson")
palette = colorRampPalette(c("steelblue","#f4f0da","darkred")) (200)
heatmap.2(x = cor.exp, col = palette,trace="none",density.info='none',dendrogram="none",Colv=FALSE,Rowv=FALSE) #for printing color key

## Beta Maturation Gene Heatmap (Fig 2C)
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) <- "sequencing2"
all_beta$sequencing2 <- factor(all_beta$sequencing2,levels=c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC', 'Fetal','Adult'))
gene.list <-c('G6PC2', 'IAPP', 'HADH', 'UCN3', 'CHGB', 'ADCYAP1', 'SIX3')
av.exp <- AverageExpression(all_beta, assays = "RNA",features  = gene.list, slot = 'data')
av.exp <- data.matrix(av.exp[[1]], rownames.force = NA)
palette = colorRampPalette(c("white","grey","black")) (10)
heatmap.2(x = av.exp,dendrogram='none',col= palette, trace="none",density.info="none",key=TRUE,scale="row",Rowv=FALSE,Colv=FALSE,cexCol=1,cexRow = 1)  

## Beta DEGs Heatmap (Fig 2D)
DefaultAssay(all_beta) = 'RNA'
all_genes = rownames(all_beta)
all_beta = ScaleData(all_beta, features = all_genes)
all_beta_shuffle = RandomSubsetData(all_beta, .999)
Idents(all_beta_shuffle) = "sequencing2"
my_levels <- c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC', 'Fetal', 'Adult')
all_beta_shuffle@active.ident <- factor(x = all_beta_shuffle@active.ident, levels = my_levels)
beta_protocol_DEGs <- FindAllMarkers(all_beta_shuffle, logfc.threshold = 0.25, only.pos = T)
beta_protocol_DEGs %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
my_levels <- c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC', 'Fetal', 'Adult')
all_beta_shuffle@active.ident <- factor(x = all_beta_shuffle@active.ident, levels = my_levels)
DoHeatmap(all_beta_shuffle, features = top20$gene, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A"))

## Downloading EC Data
setwd("")
endocrine_integrated = readRDS("endocrine_integrated.rds")
all_EC = subset(endocrine_integrated, endocrine_identity %in% "EC")
all_EC = subset(all_EC, maturation %in% c('SC-Vitro'))

## PCA & Clustering
DefaultAssay(all_EC) = "RNA"
all_EC <- NormalizeData(all_EC)
all_EC <- FindVariableFeatures(all_EC, selection.method = "vst", nfeatures = 2000)
all_EC <- ScaleData(all_EC, vars.to.regress = c('percent.mt', 'nCount_RNA', 'sequencing'))
all_EC <- RunPCA(all_EC, features = VariableFeatures(object = all_EC))
all_EC <- FindNeighbors(all_EC, dims = 1:20)
all_EC <- FindClusters(all_EC, resolution = 1)
all_EC <- RunUMAP(all_EC, dims = 1:20)

## UMAP (Fig 2E) 
DimPlot(all_EC, reduction = "umap", label = F, group.by = "sequencing",
        cols = c('#D95F02', '#7570B3', '#E7298A', '#66A61E'))+NoLegend()

## Pearson Correlation (Fig 2F)
Idents(all_EC) <-"sequencing"
DefaultAssay(all_EC) <- 'RNA'
all_EC<-FindVariableFeatures(all_EC,assay="RNA",nfeatures = 1000) 
genelist <- all_EC@assays[["RNA"]]@var.features #for most variable genes
av.exp <- AverageExpression(all_EC,features = genelist,assay="RNA",group.by ="sequencing")
av.exp <- as.matrix(av.exp[[1]])
cor.exp <- cor(av.exp,method = "pearson")
palette = colorRampPalette(c("steelblue","#f4f0da","darkred")) (200)
heatmap.2(x = cor.exp, col = palette,trace="none",density.info='none',dendrogram="none",Colv=FALSE,Rowv=FALSE) #for printing color key

## Violin Plots of EC markers (Fig 2G)
DefaultAssay(all_EC) = "RNA"
Idents(all_EC) = "sequencing"
my_levels <- c('Aug SC', 'Balboa SC', 'Veres SC', 'Weng SC')
all_EC@active.ident <- factor(x = all_EC@active.ident, levels = my_levels)
VlnPlot(all_EC, features = c("TPH1"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#66A61E'),pt.size = 0)+NoLegend()
VlnPlot(all_EC, features = c("SLC18A1"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#66A61E', '#1B9E77'),pt.size = 0)+NoLegend()
VlnPlot(all_EC, features = c("DDC"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#66A61E', '#1B9E77'),pt.size = 0)+NoLegend()
VlnPlot(all_EC, features = c("FEV"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#66A61E', '#1B9E77'),pt.size = 0)+NoLegend()
VlnPlot(all_EC, features = c("NKX6-1"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#66A61E', '#1B9E77'),pt.size = 0)+NoLegend()
VlnPlot(all_EC, features = c("LMX1A"), 
        cols = c('#D95F02', '#7570B3', '#E7298A', '#66A61E', '#1B9E77'),pt.size = 0)+NoLegend()

## EC DEGs Heatmap (Fig 2H)
DefaultAssay(all_EC) = 'RNA'
all_genes = rownames(all_EC)
all_EC = ScaleData(all_EC, features = all_genes)
all_EC_shuffle = RandomSubsetData(all_EC, .999)
Idents(all_EC_shuffle) = "sequencing"
EC_protocol_DEGs <- FindAllMarkers(all_EC_shuffle, logfc.threshold = 0.25, only.pos = T)
EC_protocol_DEGs %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
my_levels <- c('Weng SC', 'Balboa SC', 'Aug SC', 'Veres SC')
all_EC_shuffle@active.ident <- factor(x = all_EC_shuffle@active.ident, levels = my_levels)
DoHeatmap(all_EC_shuffle, features = top20$gene, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A"))
