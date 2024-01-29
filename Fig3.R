# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig 3

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

## PCA & Clustering
DefaultAssay(all_beta) = "RNA"
all_beta <- NormalizeData(all_beta)
all_beta <- FindVariableFeatures(all_beta, selection.method = "vst", nfeatures = 2000)
all_beta <- ScaleData(all_beta, vars.to.regress = c('percent.mt', 'nCount_RNA', 'sequencing', 'orig.ident'))
all_beta <- RunPCA(all_beta, features = VariableFeatures(object = all_beta))
all_beta <- FindNeighbors(all_beta, dims = 1:20)
all_beta <- FindClusters(all_beta, resolution = 1)
all_beta <- RunUMAP(all_beta, dims = 1:20)

## UMAP (Fig 3A) 
DimPlot(all_beta, reduction = "umap", label = F, group.by = "maturation",
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'))+NoLegend()

## Pearson Correlation Heatmap (Fig 3B)
Idents(all_beta) <-"maturation"
DefaultAssay(all_beta) <- 'RNA'
all_beta<-FindVariableFeatures(all_beta,assay="RNA",nfeatures = 2000) 
genelist <- all_beta@assays[["RNA"]]@var.features #for most variable genes
av.exp <- AverageExpression(all_beta,features = genelist,assay="RNA",group.by ="maturation")
av.exp <- as.matrix(av.exp[[1]])
cor.exp <- cor(av.exp,method = "pearson")
palette = colorRampPalette(c("steelblue","#f4f0da","darkred")) (20)
heatmap.2(x = cor.exp, col = palette,trace="none",density.info='none',dendrogram="none",Colv=FALSE,Rowv=FALSE) #for printing color key

## Beta DEG Heatmap (Fig 3C)
DefaultAssay(all_beta) = 'RNA'
all_genes = rownames(all_beta)
all_beta = ScaleData(all_beta, features = all_genes)
all_beta_shuffle = RandomSubsetData(all_beta, .9999)
Idents(all_beta_shuffle) = "maturation"
my_levels <- c('Adult', 'Fetal', 'SC-Vitro', 'SC-Vivo')
all_beta_shuffle@active.ident <- factor(x = all_beta_shuffle@active.ident, levels = my_levels)
adult_SC_beta_DEGs_2 <- FindAllMarkers(all_beta_shuffle, logfc.threshold = 0.25, only.pos = T)
adult_SC_beta_DEGs_2 %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50
DoHeatmap(all_beta_shuffle, features = top50$gene, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A"))

## Split beta cells for comparison
all_beta = subset(endocrine_integrated, endocrine_identity %in% "Beta")
sc_adult_beta = subset(all_beta, maturation %in% c("SC-Vitro", "Adult"))
sc_fetal_beta = subset(all_beta, maturation %in% c("SC-Vitro", "Fetal"))
adult_fetal_beta = subset(all_beta, maturation %in% c("Adult", "Fetal"))

## Find DEGs for comparison
Idents(sc_adult_beta) = "maturation"
DefaultAssay(sc_adult_beta) = "RNA"
adult_SC_beta_DEGs = FindMarkers(sc_adult_beta, ident.1 = "Adult", ident.2 = "SC-Vitro", logfc.threshold = 0)
adult_SC_beta_DEGs = subset(adult_SC_beta_DEGs, p_val_adj<0.05)

Idents(sc_fetal_beta) = "maturation"
DefaultAssay(sc_fetal_beta) = "RNA"
fetal_SC_beta_DEGs = FindMarkers(sc_fetal_beta, ident.1 = "Fetal", ident.2 = "SC-Vitro", logfc.threshold = 0)
fetal_SC_beta_DEGs = subset(fetal_SC_beta_DEGs, p_val_adj<0.05)

Idents(adult_fetal_beta) = "maturation"
DefaultAssay(adult_fetal_beta) = "RNA"
adult_fetal_beta_DEGs = FindMarkers(adult_fetal_beta, ident.1 = "Adult", ident.2 = "Fetal", logfc.threshold = 0)
adult_fetal_beta_DEGs = subset(adult_fetal_beta_DEGs, p_val_adj<0.05)

## SC-Adult Volcano (Fig 3D)
selectlabel <- c('CHGB', 'MIF', 'CPE', 'SCG5', 'SCG2', 'SCG3', 'UCN3', 'SNAP25', 'VGF', 'PTPRN', 'NUCB2', 'PDIA3', 'BMP5',
                 'NEFM', 'CRYBA2', 'C1QL1', 'CHGA', 'DCX', 'DLK1', 'LMO2', 'YWHAZ', 'SLC17A6', 'CDC42')
x <- rownames(adult_SC_beta_DEGs) %in% selectlabel
adult_SC_beta_DEGs <- rbind(adult_SC_beta_DEGs[!x,], adult_SC_beta_DEGs[x,])
keyvals <- ifelse(rownames(adult_SC_beta_DEGs) %in% selectlabel, 'green','grey')
names(keyvals)[keyvals == 'green'] <- 'picked'
names(keyvals)[keyvals == 'grey'] <- 'rest'
EnhancedVolcano(adult_SC_beta_DEGs, lab = rownames(adult_SC_beta_DEGs),selectLab=selectlabel,
                x = 'avg_log2FC', y = 'p_val_adj',xlim = c(-5, 5),title = 'Adult vs SC',
                pCutoff = NA,FCcutoff = NA, pointSize = 5,labSize = 0,colAlpha = 1.0,
                col=c('pink', 'darkgreen', 'blue', 'red3'),gridlines.minor = FALSE,
                gridlines.major = FALSE, legendPosition = 'top',drawConnectors = F,widthConnectors = 0,
                colConnectors = 'black',colCustom = keyvals, cutoffLineCol = "red")

## SC-Fetal Volcano (Fig 3D)
selectlabel <- c('IGF2', 'MEIS2', 'FOXP2', 'RBP4', 'MAFA', 'SST', 'GCK', 'CPE', 'FOXP1', 'INS',
                 'IAPP', 'TTR', 'CRYBA2', 'ONECUT2', 'FEV', 'SCG5', 'MIF', 'CHGA', 'HOPX', 'SCG3')
x <- rownames(fetal_SC_beta_DEGs) %in% selectlabel
fetal_SC_beta_DEGs <- rbind(fetal_SC_beta_DEGs[!x,], fetal_SC_beta_DEGs[x,])
keyvals <- ifelse(rownames(fetal_SC_beta_DEGs) %in% selectlabel, 'green','grey')
names(keyvals)[keyvals == 'green'] <- 'picked'
names(keyvals)[keyvals == 'grey'] <- 'rest'
EnhancedVolcano(fetal_SC_beta_DEGs, lab = rownames(fetal_SC_beta_DEGs),selectLab=selectlabel,
                x = 'avg_log2FC', y = 'p_val_adj',xlim = c(-5, 5),title = 'Fetal vs SC',
                pCutoff = NA,FCcutoff = NA, pointSize = 5,labSize = 0,colAlpha = 1.0,
                col=c('pink', 'darkgreen', 'blue', 'red3'),gridlines.minor = FALSE,
                gridlines.major = FALSE, legendPosition = 'top',drawConnectors = F,widthConnectors = 0,
                colConnectors = 'black',colCustom = keyvals, cutoffLineCol = "red")

## Adult-Fetal Volcano (Fig 3D)
selectlabel <- c('IAPP', 'MIF', 'RBP4', 'CHGB', 'SCG5', 'SHISAL2B', 'SCG2', 'ADCYAP1', 'VGF', 'G6PC2',
                 'FOXP2', 'SCGD', 'SOX4', 'DLK1', 'GCK', 'DCX', 'MAFB')
x <- rownames(adult_fetal_beta_DEGs) %in% selectlabel
adult_fetal_beta_DEGs <- rbind(adult_fetal_beta_DEGs[!x,], adult_fetal_beta_DEGs[x,])
keyvals <- ifelse(rownames(adult_fetal_beta_DEGs) %in% selectlabel, 'green','grey')
names(keyvals)[keyvals == 'green'] <- 'picked'
names(keyvals)[keyvals == 'grey'] <- 'rest'
EnhancedVolcano(adult_fetal_beta_DEGs, lab = rownames(adult_fetal_beta_DEGs),selectLab=selectlabel,
                x = 'avg_log2FC', y = 'p_val_adj',xlim = c(-5, 5),title = 'Adult vs Fetal',
                pCutoff = NA,FCcutoff = NA, pointSize = 5,labSize = 0,colAlpha = 1.0,
                col=c('pink', 'darkgreen', 'blue', 'red3'),gridlines.minor = FALSE,
                gridlines.major = FALSE, legendPosition = 'top',drawConnectors = F,widthConnectors = 0,
                colConnectors = 'black',colCustom = keyvals, cutoffLineCol = "red")

## VlnPlot of Beta markers (Fig 3E)
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
my_levels <- c('Adult', 'Fetal', 'SC-Vitro', 'SC-Vivo')
all_beta@active.ident <- factor(x = all_beta@active.ident, levels = my_levels)
VlnPlot(all_beta, features = c("INS"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()
VlnPlot(all_beta, features = c("IAPP"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()
VlnPlot(all_beta, features = c("G6PC2"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()
VlnPlot(all_beta, features = c("ABCC8"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()
VlnPlot(all_beta, features = c("SHISAL2B"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()
VlnPlot(all_beta, features = c("ADCYAP1"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()
VlnPlot(all_beta, features = c("HOPX"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()
VlnPlot(all_beta, features = c("CHGB"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()
VlnPlot(all_beta, features = c("MAFA"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()
VlnPlot(all_beta, features = c("HADH"),
        cols = c('dodgerblue', 'green3', 'magenta2', 'darkorchid4'), pt.size = 0)+NoLegend()

## Average Expression of Exocrine Markers in Beta cells (Fig 3F)
DefaultAssay(all_beta) = 'RNA'
all_genes = rownames(all_beta)
all_beta = ScaleData(all_beta, features = all_genes)
exo_genes_expression = AverageExpression(all_beta, features = c("CLPS", "CEL", "CPA1", "CTRB2"), group.by = "maturation", slot = "data")

## Heatmap of KEGG ribosomal genes (Fig 3G)
## (https://www.gsea-msigdb.org/gsea/msigdb/cards/KEGG_RIBOSOME)
KEGG_ribosome = read.delim('KEGG_RIBOSOME.v2022.1.Hs.txt')
KEGG_ribosome = KEGG_ribosome$KEGG_RIBOSOME[2:89]
Idents(all_beta_shuffle) = "maturation"
my_levels <- c('Adult', 'Fetal', 'SC-Vitro', 'SC-Vivo')
all_beta_shuffle@active.ident <- factor(x = all_beta_shuffle@active.ident, levels = my_levels)
DoHeatmap(all_beta_shuffle, features = KEGG_ribosome, slot = 'scale.data')+ scale_fill_gradientn(colors = c("blue", "white", "red"))

