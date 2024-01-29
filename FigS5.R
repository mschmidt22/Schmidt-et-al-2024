# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig S5

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

## Violin Plots of Enriched Regulons (Fig S5B)
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
my_levels <- c('SC-Vitro', 'SC-Vivo', 'Adult', 'Fetal')
all_beta@active.ident <- factor(x = all_beta@active.ident, levels = my_levels)
VlnPlot(all_beta, features = "FOXA2", pt.size = 0, split.by = "maturation", cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = "FOXA1", pt.size = 0, split.by = "maturation", cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = "ONECUT2", pt.size = 0, split.by = "maturation", cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = "PAX4", pt.size = 0, split.by = "maturation", cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = "PBX1", pt.size = 0, split.by = "maturation", cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = "SOX4", pt.size = 0, split.by = "maturation", cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = "FOXO1", pt.size = 0, split.by = "maturation", cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = "MEIS1", pt.size = 0, split.by = "maturation", cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()

## Violin Plots of Enriched Regulon Target Genes (Fig S5C)
## Regulons identified from SCENIC Analysis in Fig4

### Adding Regulon Data to seurat object
AUCmat <- AUCell::getAUC(regulonAUC)
sc_adult_fetal_beta[["AUC"]] <- CreateAssayObject(data = AUCmat)
DefaultAssay(sc_adult_fetal_beta) <- "AUC"
sc_adult_fetal_beta <- ScaleData(sc_adult_fetal_beta, assay = "AUC")

### Identification of top regulon target genes 

#### FOXA2 Top 15 regulons
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
FOXA2_targets = regulons[["FOXA2(+)"]]
FOXA2_avg_expression = AverageExpression(sc_adult_fetal_beta, features = FOXA2_targets)
FOXA2_avg_expression = FOXA2_avg_expression[1]
FOXA2_avg_expression = as.data.frame(FOXA2_avg_expression)
FOXA2_avg_expression <- FOXA2_avg_expression[1]
FOXA2_avg_expression<-FOXA2_avg_expression[order(FOXA2_avg_expression$RNA.SC.Vitro, decreasing=TRUE),,drop=F]
FOXA2_top_15_expressed_genes = rownames(FOXA2_avg_expression)
FOXA2_top_15_expressed_genes = FOXA2_top_15_expressed_genes[1:15]

#### FOXA1 Top 20 regulons
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
FOXA1_targets = regulons[["FOXA1(+)"]]
FOXA1_avg_expression = AverageExpression(sc_adult_fetal_beta, features = FOXA1_targets)
FOXA1_avg_expression = FOXA1_avg_expression[1]
FOXA1_avg_expression = as.data.frame(FOXA1_avg_expression)
FOXA1_avg_expression <- FOXA1_avg_expression[1]
FOXA1_avg_expression<-FOXA1_avg_expression[order(FOXA1_avg_expression$RNA.SC.Vitro,decreasing=TRUE),,drop=F]
FOXA1_top_20_expressed_genes = rownames(FOXA1_avg_expression)
FOXA1_top_20_expressed_genes = FOXA1_top_20_expressed_genes[1:20]

#### ONECUT2 Top 13 regulons
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
ONECUT2_targets = regulons[["ONECUT2(+)"]]
ONECUT2_avg_expression = AverageExpression(sc_adult_fetal_beta, features = ONECUT2_targets)
ONECUT2_avg_expression = ONECUT2_avg_expression[1]
ONECUT2_avg_expression = as.data.frame(ONECUT2_avg_expression)
ONECUT2_avg_expression <- ONECUT2_avg_expression[1]
ONECUT2_avg_expression<-ONECUT2_avg_expression[order(ONECUT2_avg_expression$RNA.SC.Vitro,decreasing=TRUE),,drop=F]
ONECUT2_top_13_expressed_genes = rownames(ONECUT2_avg_expression)
ONECUT2_top_13_expressed_genes = ONECUT2_top_13_expressed_genes[1:13]

#### PAX4 Top 20 regulons
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
PAX4_targets = regulons[["PAX4(+)"]]
PAX4_avg_expression = AverageExpression(sc_adult_fetal_beta, features = PAX4_targets)
PAX4_avg_expression = PAX4_avg_expression[1]
PAX4_avg_expression = as.data.frame(PAX4_avg_expression)
PAX4_avg_expression <- PAX4_avg_expression[1]
PAX4_avg_expression<-PAX4_avg_expression[order(PAX4_avg_expression$RNA.SC.Vitro,decreasing=TRUE),,drop=F]
PAX4_top_20_expressed_genes = rownames(PAX4_avg_expression)
PAX4_top_20_expressed_genes = PAX4_top_20_expressed_genes[1:20]

#### PBX1 Top 20 regulons
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
PBX1_targets = regulons[["PBX1(+)"]]
PBX1_avg_expression = AverageExpression(sc_adult_fetal_beta, features = PBX1_targets)
PBX1_avg_expression = PBX1_avg_expression[1]
PBX1_avg_expression = as.data.frame(PBX1_avg_expression)
PBX1_avg_expression <- PBX1_avg_expression[1]
PBX1_avg_expression<-PBX1_avg_expression[order(PBX1_avg_expression$RNA.SC.Vitro,decreasing=TRUE),,drop=F]
PBX1_top_20_expressed_genes = rownames(PBX1_avg_expression)
PBX1_top_20_expressed_genes = PBX1_top_20_expressed_genes[1:20]

#### SOX4 Top 20 regulons
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
SOX4_targets = regulons[["SOX4(+)"]]
SOX4_avg_expression = AverageExpression(sc_adult_fetal_beta, features = SOX4_targets)
SOX4_avg_expression = SOX4_avg_expression[1]
SOX4_avg_expression = as.data.frame(SOX4_avg_expression)
SOX4_avg_expression <- SOX4_avg_expression[1]
SOX4_avg_expression<-SOX4_avg_expression[order(SOX4_avg_expression$RNA.SC.Vitro,decreasing=TRUE),,drop=F]
SOX4_top_20_expressed_genes = rownames(SOX4_avg_expression)
SOX4_top_20_expressed_genes = SOX4_top_20_expressed_genes[1:20]

#### FOXO1 Top 10 regulons
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
FOXO1_targets = regulons[["FOXO1(+)"]]
FOXO1_avg_expression = AverageExpression(sc_adult_fetal_beta, features = FOXO1_targets)
FOXO1_avg_expression = FOXO1_avg_expression[1]
FOXO1_avg_expression = as.data.frame(FOXO1_avg_expression)
FOXO1_avg_expression <- FOXO1_avg_expression[3]
FOXO1_avg_expression<-FOXO1_avg_expression[order(FOXO1_avg_expression$RNA.Fetal,decreasing=TRUE),,drop=F]
FOXO1_top_10_expressed_genes = rownames(FOXO1_avg_expression)
FOXO1_top_10_expressed_genes = FOXO1_top_10_expressed_genes[1:10]

#### MEIS1 Top 10 regulons
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
MEIS1_targets = regulons[["MEIS1(+)"]]
MEIS1_avg_expression = AverageExpression(sc_adult_fetal_beta, features = MEIS1_targets)
MEIS1_avg_expression = MEIS1_avg_expression[1]
MEIS1_avg_expression = as.data.frame(MEIS1_avg_expression)
MEIS1_avg_expression <- MEIS1_avg_expression[3]
MEIS1_avg_expression<-MEIS1_avg_expression[order(MEIS1_avg_expression$RNA.Fetal,decreasing=TRUE),,drop=F]
MEIS1_top_16_expressed_genes = rownames(MEIS1_avg_expression)
MEIS1_top_16_expressed_genes = MEIS1_top_16_expressed_genes[1:16]

### Violin Plots
VlnPlot(all_beta, features = FOXA2_top_15_expressed_genes, split.by  = "maturation" ,stack = T, pt.size = 0, cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = FOXA1_top_20_expressed_genes, split.by  = "maturation" ,stack = T, pt.size = 0, cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = ONECUT2_top_13_expressed_genes, split.by  = "maturation" ,stack = T, pt.size = 0, cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = PAX4_top_20_expressed_genes, split.by  = "maturation" ,stack = T, pt.size = 0, cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = PBX1_top_20_expressed_genes, split.by  = "maturation" ,stack = T, pt.size = 0, cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = SOX4_top_20_expressed_genes, split.by  = "maturation" ,stack = T, pt.size = 0, cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = FOXO1_top_10_expressed_genes, split.by  = "maturation" ,stack = T, pt.size = 0, cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
VlnPlot(all_beta, features = MEIS1_top_16_expressed_genes, split.by  = "maturation" ,stack = T, pt.size = 0, cols = c('dodgerblue','green3', 'magenta2', 'darkorchid4'))+NoLegend()
