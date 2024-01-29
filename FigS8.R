# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig S8

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

# Downloading iPSC stage6wk2 SC-islet raw data from Augsornworawat et al. 2020 (DOI:10.1016/j.celrep.2020.108067)
# GSE139535 (GSM: GSM4143264)

## Download Raw Data
setwd("")
augs_data_2 = Read10X(data.dir = "")

## Creating Seurat Object
augs_2 = CreateSeuratObject(augs_data_2, project = "Augs IPSC", min.cells = 3, min.features = 200)

## Quality Control
augs_2[["percent.mt"]] = PercentageFeatureSet(augs_2, pattern = "^MT-")
VlnPlot(augs_2, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)
augs_2_QC = subset(augs_2, subset = percent.mt < 20 & nCount_RNA > 3000)

## PCA, Clustering, & Visualizing
augs_2_QC <- NormalizeData(augs_2_QC, normalization.method = "LogNormalize", scale.factor = 10000)
augs_2_QC <- FindVariableFeatures(augs_2_QC, selection.method = "vst", nfeatures = 2000)
augs_2_QC <- ScaleData(augs_2_QC)
augs_2_QC <- RunPCA(augs_2_QC, features = VariableFeatures(object = augs_2_QC))
augs_2_QC <- FindNeighbors(augs_2_QC, dims = 1:20)
augs_2_QC <- FindClusters(augs_2_QC, resolution = 1)
augs_2_QC <- RunUMAP(augs_2_QC, dims = 1:20)
DimPlot(augs_2_QC, reduction = "umap", label = T)
FeaturePlot(augs_2_QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Subsetting Endocrine Cells
augs_IPSC_endocrine = subset(augs_2_QC, idents = c(0,1,2,3,4,5,6,7,10,12,13,14))


# Downloading Hues8 stage6wk2 SC-islet single-nuclei sequencing raw data from Augsornworawat et al. 2023 (DOI: 10.1038/s41556-023-01150-8)
# GSE199636 (GSM: GSM5979674)

## Download Raw Data
augs_sn_data = Read10X_h5('GSM5979673_week2.filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)
augs_sn_data = augs_sn_data[1]
augs_sn <- CreateSeuratObject(counts = augs_sn_data$`Gene Expression`, assay = "RNA")

## Quality Control
VlnPlot(augs_sn, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)
augs_sn_qc = subset(augs_sn, subset = nCount_RNA > 1000 & nCount_RNA < 40000)

## PCA, Clustering, & Visualizing 
augs_sn_qc <- NormalizeData(augs_sn_qc, normalization.method = "LogNormalize", scale.factor = 10000)
augs_sn_qc <- FindVariableFeatures(augs_sn_qc, selection.method = "vst", nfeatures = 2000)
augs_sn_qc <- ScaleData(augs_sn_qc)
augs_sn_qc <- RunPCA(augs_sn_qc, features = VariableFeatures(object = augs_sn_qc))
augs_sn_qc <- FindNeighbors(augs_sn_qc, dims = 1:20)
augs_sn_qc <- FindClusters(augs_sn_qc, resolution = 0.5)
augs_sn_qc <- RunUMAP(augs_sn_qc, dims = 1:20)
DimPlot(augs_sn_qc, reduction = "umap", label = T)
FeaturePlot(augs_sn_qc, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Subsetting Endocrine Cells
augs_sn_endocrine = subset(augs_sn_qc, idents = c(0,1,2,3,4,5,6,7,8,910,12,13,14))


# Downloading H1 stage6wk2 SC-islet raw data from Zhu et al. 2023 (DOI: 10.1016/j.devcel.2023.03.011)
# GSE202497 (GSM: GSM6123275)

## Download Raw Data
setwd("")
zhu = readRDS('Zhu_D32.rds')

## Quality Control
VlnPlot(zhu, features = c("percent_mito","nFeature_RNA", "nCount_RNA"), ncol = 3)
zhu_qc = zhu

## PCA, Clustering, & Visualizing
zhu_qc <- NormalizeData(zhu_qc, normalization.method = "LogNormalize", scale.factor = 10000)
zhu_qc <- FindVariableFeatures(zhu_qc, selection.method = "vst", nfeatures = 2000)
zhu_qc <- ScaleData(zhu_qc)
zhu_qc <- RunPCA(zhu_qc, features = VariableFeatures(object = zhu_qc))
zhu_qc <- FindNeighbors(zhu_qc, dims = 1:20)
zhu_qc <- FindClusters(zhu_qc, resolution = 0.5)
zhu_qc <- RunUMAP(zhu_qc, dims = 1:20)
DimPlot(zhu_qc, reduction = "umap", label = T)
FeaturePlot(zhu_qc, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Subsetting Endocrine Cells
zhu_endocrine = subset(zhu_qc, idents = c(0,1,2,3,4,5,6))


# Integration of new islet datasets with adult pancreatic islets

## Merging All CHGA+ Populations
all_CHGA = merge(augs_IPSC_endocrine, c(augs_sn_endocrine, 
                                        zhu_endocrine, baron_endocrine, 
                                        fang_endocrine, xin_endocrine))

## Addition of Sequencing Column 
added_column = c()
length = ncol(all_CHGA)
count = 1
while (count <= length) {
  if(all_CHGA$orig.ident[count] %in% c("Augs IPSC")){
    value = "IPSC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("SeuratProject")){
    value = "single-nuclei"
  }
  else if(all_CHGA$orig.ident[count] %in% c("D32")){
    value = "Zhu"
  }
  else if(all_CHGA$orig.ident[count] %in% c('Adult-B1', 'Adult-B2', 'Adult-B3')){
    value = "Baron"
  }
  else if(all_CHGA$orig.ident[count] %in% c('Adult-F1', 'Adult-F2', 'Adult-F3', 'Adult-F4', 'Adult-F6')){
    value = "Fang"
  }
  else{
    value = "Xin"
  }
  added_column = append(added_column, value)
  count = count + 1
}

View(added_column)
all_CHGA$sequencing = added_column

## Addition of Maturation Column 
added_column = c()
length = ncol(all_CHGA)
count = 1
while (count <= length) {
  if(all_CHGA$orig.ident[count] %in% c("Augs IPSC")){
    value = "IPSC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("SeuratProject")){
    value = "single-nuclei"
  }
  else if(all_CHGA$orig.ident[count] %in% c("D32")){
    value = "Zhu"
  }
  else{
    value = "Adult"
  }
  added_column = append(added_column, value)
  count = count + 1
}

View(added_column)
all_CHGA$maturation = added_column

## Removing all mitochondrial genes (MT) from each dataset for consistency 
DefaultAssay(all_CHGA) = "RNA"
all_genes = rownames(all_CHGA)
anti_mt_genes = all_genes[grep(x = all_genes, pattern ="^MT-", invert = T)]
mt_genes = all_genes[grepl(x = all_genes, pattern = "MT-")]
all_CHGA = subset(all_CHGA, features = anti_mt_genes)

## Setting Up Object for Integration
split.list <- SplitObject(all_CHGA, split.by = "sequencing")
split.list <- lapply(X = split.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = split.list)

## Performing Integration
endocrine.anchors <- FindIntegrationAnchors(object.list = split.list, anchor.features = features)
CHGA_2_integrated <- IntegrateData(anchorset = endocrine.anchors)
DefaultAssay(CHGA_2_integrated) <- "integrated"
CHGA_2_integrated <- ScaleData(CHGA_2_integrated, vars.to.regress = c('nCount_RNA'))
CHGA_2_integrated <- RunPCA(CHGA_2_integrated, npcs = 30, verbose = FALSE)
CHGA_2_integrated <- RunUMAP(CHGA_2_integrated, reduction = "pca", dims = 1:30)
CHGA_2_integrated <- FindNeighbors(CHGA_2_integrated, reduction = "pca", dims = 1:30)
CHGA_2_integrated <- FindClusters(CHGA_2_integrated, resolution = 2)

## Cell Typing
CHGA_2_integrated = RenameIdents(CHGA_2_integrated, '0' = "EC", '1' = "Alpha",
                                 '2' = "EC", '3' = "Beta", '4' = "Alpha", 
                                 '5' = "Delta", '6'="Beta", '7'="Alpha",
                                 '8'="Beta",'9'="Beta", '10'="EC", '11'="Beta",
                                 '12'="Alpha",'13' = "Beta", '14'="Alpha", '15'="Alpha",
                                 '16'="Alpha",'17'= "Alpha", '18' = "Beta", '19' = "PP", '20' = "Alpha"
                                 , '21' = "Delta", '22' = "Alpha", '23' = "EC", '24' = "Beta"
                                 , '25' = "Beta", '26' = "Alpha", '27' = "Alpha", '28' = "Alpha"
                                 , '29' = "EC", '30' = "Delta", '31' = "Beta", '32' = "Exo", '33' = "Beta"
                                 , '34' = "Epsilon")
CHGA_2_integrated$CHGA_identity = CHGA_2_integrated@active.ident

## UMAPs (Fig S8A & S8B)
DimPlot(CHGA_2_integrated, reduction = "umap", repel = TRUE, label = F, cols = c("#4DAF4A", "#377EB8", "#E41A1C", "#984EA3",
                                                                                 "#7570B3", "#A6761D", "#E7298A"))+NoLegend()
DimPlot(CHGA_2_integrated, reduction = "umap", repel = TRUE, label = T, split.by = "maturation")

## Neural Identity Violin Panels (Fig S8C)
beta = subset(CHGA_2_integrated, CHGA_identity %in% "Beta")
DefaultAssay(beta) = "RNA"
Idents(beta) = "maturation"
VlnPlot(beta, features = c('KIF5C', 'BASP1', 'CD24', 'MARCKS', 'MAP1B', 'MDK',
                           'MLLT11', 'RUFY3', 'NEFM', 'HMGCR', 'NAV1', 'NTM', 'PLPPR5', 'CHODL', 'EVL', 
                           'AKT1', 'ABLIM1', 'GAP43', 'NEFL'), 
        split.by  = "maturation",stack = T, 
        pt.size = 0, cols = c('#C77CFF', '#F8766D', '#7CAE00', '#00BFC4'))+NoLegend()
VlnPlot(beta, features = c('H3F3A', 'ISL1', 'GTF2I', 'ONECUT2', 'DPYSL2', 'ASCL2', 'SOX4', 'SOX11', 'RUNX1T1', 'ARID1A',
                           'ARID1B', 'PBX1', 'ELAVL2', 'ASCL1'), 
        split.by  = "maturation",stack = T, 
        pt.size = 0, cols = c('#C77CFF', '#F8766D', '#7CAE00', '#00BFC4'))+NoLegend()
VlnPlot(beta, features = c('STMN1', 'STMN2', 'DCX', 'PLXNA2','VLDLR', 'CRMP1', 'LSAMP', 
                           'NCAM1', 'PLXNC1', 'CXCL12','STMN4'), 
        split.by  = "maturation", stack = T, 
        pt.size = 0, cols = c('#C77CFF', '#F8766D', '#7CAE00', '#00BFC4'))+NoLegend()
VlnPlot(beta, features = c('MARCKSL1', 'C1QL1', 'RAB3B', 'CADPS', 'CAMK2N1', 'NBEA', 'SYT11',
                           'APOE', 'CALY','CHRNA3'), 
        split.by  = "maturation", stack = T, 
        pt.size = 0, cols = c('#C77CFF', '#F8766D', '#7CAE00', '#00BFC4'))+NoLegend()
VlnPlot(beta, features = c('CALB2', 'FXYD6', 'PCP4', 'ATP2B1', 'CELF4', 'CACNA2D1', 'CACNA1A'), 
        split.by  = "maturation",stack = T, 
        pt.size = 0, cols = c('#C77CFF', '#F8766D', '#7CAE00', '#00BFC4'))+NoLegend()
