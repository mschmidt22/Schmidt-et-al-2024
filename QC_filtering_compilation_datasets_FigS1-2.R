# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# This code was used to download published datasets, align gene names, quality control, & perform PCA/Clustering/Cell typing of individual datasets
# Fig S1 & Fig S2

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

# Downloading Hues8 stage6wk2 SC-islet raw data from Augsornworawat et al. 2020 (DOI:10.1016/j.celrep.2020.108067)
# GSE151117 (GSM: GSM4567006)

## Download Raw Data
setwd("")
augs_data = Read10X(data.dir = "")

## Editing Gene Names
rownames(augs_data)[rownames(augs_data) == "ATP5I"] ="ATP5ME"
rownames(augs_data)[colnames(augs_data) == "ATP5G2"] ="ATP5MC2"
rownames(augs_data)[rownames(augs_data) == "ATP5E"] ="ATP5F1E"
rownames(augs_data)[rownames(augs_data) == "ATP5L"] ="ATP5MG"
rownames(augs_data)[rownames(augs_data) == "TCEB2"] ="ELOB"
rownames(augs_data)[rownames(augs_data) == "SHFM1"] ="SEM1"
rownames(augs_data)[rownames(augs_data) == "USMG5"] ="ATP5MD"
rownames(augs_data)[rownames(augs_data) == "NGFRAP1"] ="BEX3"
rownames(augs_data)[rownames(augs_data) == "FAM213A"] ="PRXL2A"
rownames(augs_data)[rownames(augs_data) == "GNB2L1"] ="RACK1"
rownames(augs_data)[rownames(augs_data) == "TCEB3"] ="ELOA"
rownames(augs_data)[rownames(augs_data) == "C8orf59"] ="RBIS"
rownames(augs_data)[rownames(augs_data) == "C19orf24"] ="FAM174C"
rownames(augs_data)[rownames(augs_data) == "C20orf24"] ="RAB5IF"
rownames(augs_data)[rownames(augs_data) == "FAM173A"] ="ANTKMT"
rownames(augs_data)[rownames(augs_data) == "C6orf1"] ="SMIM29"
rownames(augs_data)[rownames(augs_data) == "SEPT2"] ="SEPTIN2"
rownames(augs_data)[rownames(augs_data) == "FAM57B"] ="TLCD3B"
rownames(augs_data)[rownames(augs_data) == "PVRL2"] ="NECTIN2"
rownames(augs_data)[rownames(augs_data) == "FAM96A"] ="CIAO2A"
rownames(augs_data)[rownames(augs_data) == "KIAA1715"] ="LNPK"
rownames(augs_data)[rownames(augs_data) == "RTFDC1"] ="RTF2"
rownames(augs_data)[rownames(augs_data) == "FAM134A"] ="RETREG2"
rownames(augs_data)[rownames(augs_data) == "STRA13"] ="CENPX"
rownames(augs_data)[rownames(augs_data) == "C9orf3"] ="AOPEP"
rownames(augs_data)[rownames(augs_data) == "HBT8"] ="PWAR6"
rownames(augs_data)[rownames(augs_data) == "PVRL3"] ="NECTIN3"
rownames(augs_data)[rownames(augs_data) == "WRB"] ="GET1"
rownames(augs_data)[rownames(augs_data) == "ASNA1"] ="GET3"
rownames(augs_data)[rownames(augs_data) == "FAM60A"] ="SINHCAF"
rownames(augs_data)[rownames(augs_data) == "FAM63B"] ="MINDY2"
rownames(augs_data)[rownames(augs_data) == "COL4A3BP"] ="CERT1"
rownames(augs_data)[rownames(augs_data) == "WHSC1"] ="NSD2"
rownames(augs_data)[rownames(augs_data) == "APOPT1"] ="COA8"
rownames(augs_data)[rownames(augs_data) == "TOMM70A"] ="TOMM70"
rownames(augs_data)[rownames(augs_data) == "PVRL1"] ="NECTIN1"
rownames(augs_data)[rownames(augs_data) == "TMEM261"] ="DMAC1"
rownames(augs_data)[rownames(augs_data) == "FAM195A"] ="MCRIP2"
rownames(augs_data)[rownames(augs_data) == "PPP2R4"] ="PTPA"
rownames(augs_data)[rownames(augs_data) == "SUMO2.1"] ="SUMO2"
rownames(augs_data)[rownames(augs_data) == "UQR11.1"] ="UQCR11"
rownames(augs_data)[rownames(augs_data) == "C11orf73"] ="HIKESHI"
rownames(augs_data)[rownames(augs_data) == "C19orf60"] ="REX1BD"
rownames(augs_data)[rownames(augs_data) == "ZCCHC11"] ="TUT4"
rownames(augs_data)[rownames(augs_data) == "LHFP"] ="LHFPL6"
rownames(augs_data)[rownames(augs_data) == "MESDC2"] ="MESD"
rownames(augs_data)[rownames(augs_data) == "FAM127B"] ="RTL8A"
rownames(augs_data)[rownames(augs_data) == "C7orf55"] ="FMC1"
rownames(augs_data)[rownames(augs_data) == "LINC00493"] ="SMIM26"
rownames(augs_data)[rownames(augs_data) == "WBSCR22"] ="BUD23"
rownames(augs_data)[rownames(augs_data) == "FAM105A"] ="OTULINL"
rownames(augs_data)[rownames(augs_data) == "FAM127A"] ="RTL8C"
rownames(augs_data)[rownames(augs_data) == "C14orf1"] ="ERG28"
rownames(augs_data)[rownames(augs_data) == "SEPT2"] ="SEPTIN2"
rownames(augs_data)[rownames(augs_data) == "VIMP"] ="SELENOS"
rownames(augs_data)[rownames(augs_data) == "C16orf13"] ="METTL26"
rownames(augs_data)[rownames(augs_data) == "SELM"] ="SELENOM"
rownames(augs_data)[rownames(augs_data) == "FAM195B"] ="MCRIP1"
rownames(augs_data)[rownames(augs_data) == "MLLT4"] ="AFDN"
rownames(augs_data)[rownames(augs_data) == "WBP5"] ="TCEAL9"
rownames(augs_data)[rownames(augs_data) == "SEPT11"] ="SEPTIN11"
rownames(augs_data)[rownames(augs_data) == "C19orf70"] ="MICOS13"
rownames(augs_data)[rownames(augs_data) == "WHSC1L1"] ="NSD3"
rownames(augs_data)[rownames(augs_data) == "ATP5H"] ="ATP5PD"
rownames(augs_data)[rownames(augs_data) == "TCEB1"] ="ELOC"
rownames(augs_data)[rownames(augs_data) == "ATP5C1"] ="ATP5F1C"
rownames(augs_data)[rownames(augs_data) == "SEP15"] ="SELENOF"
rownames(augs_data)[rownames(augs_data) == "FAM96B"] ="CIA02B"
rownames(augs_data)[rownames(augs_data) == "C17orf89"] ="NDUFAF8"
rownames(augs_data)[rownames(augs_data) == "APOA1BP"] ="NAXE"
rownames(augs_data)[rownames(augs_data) == "ATP5F1"] ="ATP5PB"
rownames(augs_data)[rownames(augs_data) == "FAM46A"] ="TENT5A"
rownames(augs_data)[rownames(augs_data) == "C14orf166"] ="RTRAF"
rownames(augs_data)[rownames(augs_data) == "SEPT7"] ="SEPTIN7"
rownames(augs_data)[rownames(augs_data) == "MINOS1"] ="MICOS10"
rownames(augs_data)[rownames(augs_data) == "C11orf31"] ="SELENOH"
rownames(augs_data)[rownames(augs_data) == "GLTSCR2"] ="NOP53"
rownames(augs_data)[rownames(augs_data) == "LINC01420"] ="NBDY"
rownames(augs_data)[rownames(augs_data) == "FAM159B"] ="SHISAL2B"
rownames(augs_data)[rownames(augs_data) == "C7orf73"] ="STMP1"
rownames(augs_data)[rownames(augs_data) == "C19orf43"] ="TRIR"
rownames(augs_data)[rownames(augs_data) == "SELT"] ="SELENOT"
rownames(augs_data)[rownames(augs_data) == "AES"] ="TLE5"
rownames(augs_data)[rownames(augs_data) == "ATP5G1"] ="ATP5MC1"
rownames(augs_data)[rownames(augs_data) == "ATP5B"] ="ATP5F1B"
rownames(augs_data)[rownames(augs_data) == "SEPW1"] ="SELENOW"
rownames(augs_data)[rownames(augs_data) == "ATP5D"] ="ATP5F1D"
rownames(augs_data)[rownames(augs_data) == "ATP5O"] ="ATP5PO"
rownames(augs_data)[rownames(augs_data) == "SELK"] ="SELENOK"
rownames(augs_data)[rownames(augs_data) == "HN1"] ="JPT1"
rownames(augs_data)[rownames(augs_data) == "ATPIF1"] ="ATP5IF1"
rownames(augs_data)[rownames(augs_data) == "ATP5J"] ="ATP5PF"
rownames(augs_data)[rownames(augs_data) == "ATP5A1"] ="ATP5F1A"
rownames(augs_data)[rownames(augs_data) == "ATP5G3"] ="ATP5MC3"
rownames(augs_data)[rownames(augs_data) == "ATP5J2"] ="ATP5MF"
rownames(augs_data)[rownames(augs_data) == "MYEOV2"] ="COPS9"
rownames(augs_data)[rownames(augs_data) == "UQCR11.1"] ="UQCR11"
rownames(augs_data)[rownames(augs_data) == "SUMO2.1"] ="SUMO2"
rownames(augs_data)[rownames(augs_data) == "LRRC75A-AS1"] ="SNHG29"
rownames(augs_data)[rownames(augs_data) == "SEPT4"] ="SEPTIN4"
rownames(augs_data)[rownames(augs_data) == "LARGE1"] ="LARGE"

## Creating Seurat Object
augs = CreateSeuratObject(augs_data, project = "Augsornworawat", min.cells = 3, min.features = 200)

## Quality Control
augs[["percent.mt"]] = PercentageFeatureSet(augs, pattern = "^MT-")
VlnPlot(augs, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                        axis.title.y = element_text(size=0),
                                                                                        axis.text.x = element_text(size=0),
                                                                                        axis.text.y = element_text(size=20),
                                                                                        title = element_text(size=12))
augs_QC = subset(augs, subset = percent.mt < 25 & nFeature_RNA > 2500 & nCount_RNA < 40000 & nCount_RNA > 8000)

## Adding Metadata
original_dataset = rep("Aug_SC_vitro_sample_1", 8615)
age = rep("SC-Vitro", 8615)
BMI = rep("SC-Vitro", 8615)
gender = rep(NA, 8615)
augs_QC$original_dataset = original_dataset
augs_QC$age = age
augs_QC$BMI = BMI
augs_QC$gender = gender

## PCA & Clustering
augs_QC <- NormalizeData(augs_QC, normalization.method = "LogNormalize", scale.factor = 10000)
augs_QC <- FindVariableFeatures(augs_QC, selection.method = "vst", nfeatures = 2000)
augs_QC <- ScaleData(augs_QC)
augs_QC <- RunPCA(augs_QC, features = VariableFeatures(object = augs_QC))
augs_QC <- FindNeighbors(augs_QC, dims = 1:20)
augs_QC <- FindClusters(augs_QC, resolution = 1)
augs_QC <- RunUMAP(augs_QC, dims = 1:20)
DimPlot(augs_QC, reduction = "umap", label = T)
FeaturePlot(augs_QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
augs_QC = RenameIdents(augs_QC, '0' = "Alpha", '1' = "Alpha",
                       '2' = "Alpha", '3' = "Beta", '4' = "Alpha", 
                       '5' = "Alpha", '6'="Alpha", '7'="Alpha",
                       '8'="EC",'9'="EC", '10'="Prolif", '11' = "Alpha",
                       '12' = "Alpha", '13' = "Prolif")
augs_QC$first_identity = augs_QC@active.ident

## UMAP
DimPlot(augs_QC, reduction = "umap", label = F, cols = c("#377EB8", "#E41A1C", "#4DAF4A", "#7570B3"))+NoLegend()

## Subsetting CHGA+ Cells
augs_endocrine = subset(augs_QC, identity = c("Beta", "Alpha", "EC", "Prolif"))


# Downloading Hues8 stage6wk2 SC-islet raw data from Veres et al. 2019 (DOI:10.1038/s41586-019-1168-5)
# GSE114412 (GSM: GSM3141964, GSM3141970, & GSM3141976)

## Download Raw Data
setwd("")
veres_raw<-read.table("",header=T,sep="\t")
rownames(veres_raw)<- veres_raw[,1]
veres_raw <- veres_raw[,-1]

## Editing Gene Names
rownames(veres_raw)[rownames(veres_raw) == "ATP5I"] ="ATP5ME"
rownames(veres_raw)[rownames(veres_raw) == "ATP5G2"] ="ATP5MC2"
rownames(veres_raw)[rownames(veres_raw) == "ATP5E"] ="ATP5F1E"
rownames(veres_raw)[rownames(veres_raw) == "ATP5L"] ="ATP5MG"
rownames(veres_raw)[rownames(veres_raw) == "USMG5"] ="ATP5MD"
rownames(veres_raw)[rownames(veres_raw) == "FAM213A"] ="PRXL2A"
rownames(veres_raw)[rownames(veres_raw) == "C8orf59"] ="RBIS"
rownames(veres_raw)[rownames(veres_raw) == "C19orf24"] ="FAM174C"
rownames(veres_raw)[rownames(veres_raw) == "C20orf24"] ="RAB5IF"
rownames(veres_raw)[rownames(veres_raw) == "FAM173A"] ="ANTKMT"
rownames(veres_raw)[rownames(veres_raw) == "C6orf1"] ="SMIM29"
rownames(veres_raw)[rownames(veres_raw) == "SEPT2"] ="SEPTIN2"
rownames(veres_raw)[rownames(veres_raw) == "FAM57B"] ="TLCD3B"
rownames(veres_raw)[rownames(veres_raw) == "FAM96A"] ="CIAO2A"
rownames(veres_raw)[rownames(veres_raw) == "RTFDC1"] ="RTF2"
rownames(veres_raw)[rownames(veres_raw) == "FAM134A"] ="RETREG2"
rownames(veres_raw)[rownames(veres_raw) == "C9orf3"] ="AOPEP"
rownames(veres_raw)[rownames(veres_raw) == "WRB"] ="GET1"
rownames(veres_raw)[rownames(veres_raw) == "ASNA1"] ="GET3"
rownames(veres_raw)[rownames(veres_raw) == "FAM60A"] ="SINHCAF"
rownames(veres_raw)[rownames(veres_raw) == "COL4A3BP"] ="CERT1"
rownames(veres_raw)[rownames(veres_raw) == "APOPT1"] ="COA8"
rownames(veres_raw)[rownames(veres_raw) == "TMEM261"] ="DMAC1"
rownames(veres_raw)[rownames(veres_raw) == "C19orf60"] ="REX1BD"
rownames(veres_raw)[rownames(veres_raw) == "ZCCHC11"] ="TUT4"
rownames(veres_raw)[rownames(veres_raw) == "LHFP"] ="LHFPL6"
rownames(veres_raw)[rownames(veres_raw) == "MESDC2"] ="MESD"
rownames(veres_raw)[rownames(veres_raw) == "FAM127B"] ="RTL8A"
rownames(veres_raw)[rownames(veres_raw) == "WBSCR22"] ="BUD23"
rownames(veres_raw)[rownames(veres_raw) == "FAM105A"] ="OTULINL"
rownames(veres_raw)[rownames(veres_raw) == "FAM127A"] ="RTL8C"
rownames(veres_raw)[rownames(veres_raw) == "C14orf1"] ="ERG28"
rownames(veres_raw)[rownames(veres_raw) == "SEPT2"] ="SEPTIN2"
rownames(veres_raw)[rownames(veres_raw) == "SEPT11"] ="SEPTIN11"
rownames(veres_raw)[rownames(veres_raw) == "C19orf70"] ="MICOS13"
rownames(veres_raw)[rownames(veres_raw) == "ATP5H"] ="ATP5PD"
rownames(veres_raw)[rownames(veres_raw) == "ATP5C1"] ="ATP5F1C"
rownames(veres_raw)[rownames(veres_raw) == "FAM96B"] ="CIA02B"
rownames(veres_raw)[rownames(veres_raw) == "ATP5F1"] ="ATP5PB"
rownames(veres_raw)[rownames(veres_raw) == "FAM46A"] ="TENT5A"
rownames(veres_raw)[rownames(veres_raw) == "C14orf166"] ="RTRAF"
rownames(veres_raw)[rownames(veres_raw) == "SEPT7"] ="SEPTIN7"
rownames(veres_raw)[rownames(veres_raw) == "MINOS1"] ="MICOS10"
rownames(veres_raw)[rownames(veres_raw) == "GLTSCR2"] ="NOP53"
rownames(veres_raw)[rownames(veres_raw) == "FAM159B"] ="SHISAL2B"
rownames(veres_raw)[rownames(veres_raw) == "C7orf73"] ="STMP1"
rownames(veres_raw)[rownames(veres_raw) == "C19orf43"] ="TRIR"
rownames(veres_raw)[rownames(veres_raw) == "AES"] ="TLE5"
rownames(veres_raw)[rownames(veres_raw) == "ATP5G1"] ="ATP5MC1"
rownames(veres_raw)[rownames(veres_raw) == "ATP5B"] ="ATP5F1B"
rownames(veres_raw)[rownames(veres_raw) == "ATP5D"] ="ATP5F1D"
rownames(veres_raw)[rownames(veres_raw) == "ATP5O"] ="ATP5PO"
rownames(veres_raw)[rownames(veres_raw) == "HN1"] ="JPT1"
rownames(veres_raw)[rownames(veres_raw) == "ATPIF1"] ="ATP5IF1"
rownames(veres_raw)[rownames(veres_raw) == "ATP5J"] ="ATP5PF"
rownames(veres_raw)[rownames(veres_raw) == "ATP5A1"] ="ATP5F1A"
rownames(veres_raw)[rownames(veres_raw) == "ATP5G3"] ="ATP5MC3"
rownames(veres_raw)[rownames(veres_raw) == "ATP5J2"] ="ATP5MF"
rownames(veres_raw)[rownames(veres_raw) == "C14orf2"] ="ATP5MPL"
rownames(veres_raw)[rownames(veres_raw) == "SEPT4"] ="SEPTIN4"
rownames(veres_raw)[rownames(veres_raw) == "LARGE1"] ="LARGE"

## Creating Seurat Object
veres <-CreateSeuratObject(counts=veres_raw, project="veres",min.cells = 0, min.features = 20)

## Quality Control
veres[["percent.mt"]] = PercentageFeatureSet(veres, pattern = "^MT-")
VlnPlot(veres, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                         axis.title.y = element_text(size=0),
                                                                                         axis.text.x = element_text(size=0),
                                                                                         axis.text.y = element_text(size=20),
                                                                                         title = element_text(size=12))
veres_QC = subset(veres, subset = nFeature_RNA > 300 & nCount_RNA < 7500 & nCount_RNA > 1500)

## Adding Metadata
original_dataset = rep("Veres_SC_vitro_sample_1", 4577)
age = rep("SC-Vitro", 4577)
BMI = rep("SC-Vitro", 4577)
gender = rep(NA, 4577)
veres_QC$original_dataset = original_dataset
veres_QC$age = age
veres_QC$BMI = BMI
veres_QC$gender = gender

## PCA & Clustering
veres_QC <- NormalizeData(veres_QC, normalization.method = "LogNormalize", scale.factor = 10000)
veres_QC <- FindVariableFeatures(veres_QC, selection.method = "vst", nfeatures = 2000)
veres_QC <- ScaleData(veres_QC)
veres_QC <- RunPCA(veres_QC, features = VariableFeatures(object = veres_QC))
veres_QC <- FindNeighbors(veres_QC, dims = 1:20)
veres_QC <- FindClusters(veres_QC, resolution = 2.5)
veres_QC <- RunUMAP(veres_QC, dims = 1:20)
DimPlot(veres_QC, reduction = "umap", label = T)
FeaturePlot(veres_QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
veres_QC = RenameIdents(veres_QC, '0' = "Beta", '1' = "Beta",
                        '2' = "EC", '3' = "Alpha", '4' = "EC", 
                        '5' = "Beta", '6'="EC", '7'="Beta",
                        '8'="Alpha",'9'="CHGA-", '10'="Beta", '11' = "Alpha",
                        '12' = "EC", '13' = "Beta", '14' = "EC", '15' = "Alpha",
                        '16' = "Exo (CHGA+)", '17' = "CHGA-", '18' = "EC", '19' = "CHGA-", '20' = "Delta",
                        '21' = "Prolif", '22' = "CHGA-")
veres_QC$first_identity = veres_QC@active.ident

## UMAP
DimPlot(veres_QC, reduction = "umap", label = F, cols = c("#E41A1C", "#4DAF4A", "#377EB8", "grey", "yellow2", "#984EA3", "#7570B3"))+NoLegend()

## Subsetting CHGA+ Cells
veres_endocrine = subset(veres_QC, identity %in% c('Beta', 'Alpha', 'Delta', 'EC', 'Prolif', 'Exo (CHGA+)'))


# Downloading H1 stage7 SC-islet raw data from Weng et al. 2020 (DOI:10.1038/s42255-020-00314-2)
# GSE143783 (GSM: GSM4274582)

## Download Raw Data
setwd("")
weng_raw<-read.table("",header=T,sep="\t")
rownames(weng_raw)<- weng_raw[,1]
weng_raw <- weng_raw[,-1]

## Editing Gene Names
rownames(weng_raw)[rownames(weng_raw) == "H3-3B"] ="H3F3B"
rownames(weng_raw)[rownames(weng_raw) == "H2AJ"] ="H2AFJ"
rownames(weng_raw)[rownames(weng_raw) == "H2AZ1"] ="H2AFZ"
rownames(weng_raw)[rownames(weng_raw) == "MACROH2A1"] ="H2AFY"
rownames(weng_raw)[rownames(weng_raw) == "H1-0"] ="H1F0"
rownames(weng_raw)[rownames(weng_raw) == "EPRS1"] ="EPRS"
rownames(weng_raw)[rownames(weng_raw) == "PME3IP1"] ="FAM192A"
rownames(weng_raw)[rownames(weng_raw) == "TARS1"] ="TARS"
rownames(weng_raw)[rownames(weng_raw) == "H1-10"] ="H1FX"
rownames(weng_raw)[rownames(weng_raw) == "H4C3"] ="HIST1H4C"
rownames(weng_raw)[rownames(weng_raw) == "NARS1"] ="NARS"
rownames(weng_raw)[rownames(weng_raw) == "H2AZ2"] ="H2AFV"
rownames(weng_raw)[rownames(weng_raw) == "H1-4"] ="HIST1H1E"

## Creating Seurat Object
weng <-CreateSeuratObject(counts=weng_raw, project="weng",min.cells = 0, min.features = 20)

## Quality Control
weng[["percent.mt"]] = PercentageFeatureSet(weng, pattern = "^MT-")
VlnPlot(weng, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                        axis.title.y = element_text(size=0),
                                                                                        axis.text.x = element_text(size=0),
                                                                                        axis.text.y = element_text(size=20),
                                                                                        title = element_text(size=12))
weng_QC = subset(weng, subset = nFeature_RNA > 200 & percent.mt < 25 & nCount_RNA < 2000 & nCount_RNA > 200)

## Adding Metadata
original_dataset = rep("Weng_SC_vitro_sample_1", 9229)
age = rep("SC-Vitro", 9229)
BMI = rep("SC-Vitro", 9229)
gender = rep(NA, 9229)
weng_QC$original_dataset = original_dataset
weng_QC$age = age
weng_QC$BMI = BMI
weng_QC$gender = gender

## PCA & Clustering
weng_QC <- NormalizeData(weng_QC, normalization.method = "LogNormalize", scale.factor = 10000)
weng_QC <- FindVariableFeatures(weng_QC, selection.method = "vst", nfeatures = 2000)
weng_QC <- ScaleData(weng_QC)
weng_QC <- RunPCA(weng_QC, features = VariableFeatures(object = weng_QC))
weng_QC <- FindNeighbors(weng_QC, dims = 1:20)
weng_QC <- FindClusters(weng_QC, resolution = 0.7)
weng_QC <- RunUMAP(weng_QC, dims = 1:20)
DimPlot(weng_QC, reduction = "umap", label = T)
FeaturePlot(weng_QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
weng_QC = RenameIdents(weng_QC, '0' = "CHGA-", '1' = "CHGA-",
                       '2' = "CHGA-", '3' = "EC", '4' = "CHGA-", 
                       '5' = "Alpha", '6'="Beta", '7'="CHGA-",
                       '8'="CHGA-",'9'="CHGA-", '10'="CHGA-", '11' = "NeuroEndo",
                       '12' = "NeuroEndo", '13' = "Prolif")
weng_QC$first_identity = weng_QC@active.ident

## UMAP
DimPlot(weng_QC, reduction = "umap", label = F, cols = c("grey", "#4DAF4A", "#377EB8", "#E41A1C", "#A6761D", "#7570B3"))+NoLegend()

## Subsetting CHGA+ Cells
weng_endocrine = subset(weng_QC, identity = c('Beta', 'EC', 'Alpha', 'NeuroEndo', 'Prolif'))


# Downloading H1 stage7wk0 SC-islet raw data from Balboa et al. 2022 (DOI:10.1038/s41587-022-01219-z)
# GSE167880 (GSM: GSM5114461, GSM5114462, & GSM5114463)

## Download Raw Data & Creating Seurat Object
setwd("")
balboa_data_1 = Read10X(data.dir = "")
balboa_1 = CreateSeuratObject(balboa_data_1, project = "Balboa", min.cells = 3, min.features = 200)
original_dataset = rep("Balboa_SC_vitro_sample_1", 4595)
balboa_1$original_dataset = original_dataset

balboa_data_2 = Read10X(data.dir = "")
balboa_2 = CreateSeuratObject(balboa_data_2, project = "Balboa", min.cells = 3, min.features = 200)
original_dataset = rep("Balboa_SC_vitro_sample_2", 5763)
balboa_2$original_dataset = original_dataset

balboa_data_3 = Read10X(data.dir = "")
balboa_3 = CreateSeuratObject(balboa_data_3, project = "Balboa", min.cells = 3, min.features = 200)
original_dataset = rep("Balboa_SC_vitro_sample_3", 4030)
balboa_3$original_dataset = original_dataset

balboa = merge(balboa_1, c(balboa_2, balboa_3))

## Quality Control
balboa[["percent.mt"]] = PercentageFeatureSet(balboa, pattern = "^MT-")
VlnPlot(balboa, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                          axis.title.y = element_text(size=0),
                                                                                          axis.text.x = element_text(size=0),
                                                                                          axis.text.y = element_text(size=20),
                                                                                          title = element_text(size=12))
balboa_QC = subset(balboa, subset = nFeature_RNA > 2000 & percent.mt < 25 & nCount_RNA < 26000)

## Adding Metadata
age = rep("SC-Vitro", 8641)
BMI = rep("SC-Vitro", 8641)
gender = rep(NA, 8641)
balboa_QC$age = age
balboa_QC$BMI = BMI
balboa_QC$gender = gender

## PCA & Clustering
balboa_QC <- NormalizeData(balboa_QC, normalization.method = "LogNormalize", scale.factor = 10000)
balboa_QC <- FindVariableFeatures(balboa_QC, selection.method = "vst", nfeatures = 2000)
balboa_QC <- ScaleData(balboa_QC)
balboa_QC <- RunPCA(balboa_QC, features = VariableFeatures(object = balboa_QC))
balboa_QC <- FindNeighbors(balboa_QC, dims = 1:20)
balboa_QC <- FindClusters(balboa_QC, resolution = 1)
balboa_QC <- RunUMAP(balboa_QC, dims = 1:20)
DimPlot(balboa_QC, reduction = "umap", label = T)
FeaturePlot(balboa_QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
balboa_QC = RenameIdents(balboa_QC, '0' = "EC", '1' = "Beta",
                         '2' = "Alpha", '3' = "Beta", '4' = "Beta", 
                         '5' = "Beta", '6'="EC", '7'="Alpha",
                         '8'="CHGA-",'9'="EC", '10'="NeuroEndo", '11' = "CHGA-",
                         '12' = "Delta", '13' = "Prolif", '14' = "CHGA-", '15' = "Exo (CHGA+)")
balboa_QC$first_identity = balboa_QC@active.ident

## UMAP
DimPlot(balboa_QC, reduction = "umap", label = F, cols = c("#4DAF4A", "#E41A1C", "#377EB8", "grey", "#A6761D", "#984EA3", "#7570B3", "yellow2"))+NoLegend()

## Subsetting CHGA+ Cells
balboa_endocrine = subset(balboa_QC, identity = c('EC', 'Beta', 'Alpha', 'NeuroEndo', 'Delta', 'Prolif', 'Exo (CHGA+)'))


# Downloading H1 1-month TXP SC-islet raw data from Balboa et al. 2022 (DOI:10.1038/s41587-022-01219-z)
# GSE167880 (GSM: GSM5114469, GSM5114470, & GSM5114471)

## Download Raw Data & Creating Seurat Object
setwd("")
balboa_vivo_1_data_1 = Read10X(data.dir = "")
balboa_vivo_1_1 = CreateSeuratObject(balboa_vivo_1_data_1, project = "Balboa vivo 1 month", min.cells = 3, min.features = 200)
original_dataset = rep("Balboa_SC_vivo_1_month_sample_1", 1772)
balboa_vivo_1_1$original_dataset = original_dataset

balboa_vivo_1_data_2 = Read10X(data.dir = "")
balboa_vivo_1_2 = CreateSeuratObject(balboa_vivo_1_data_2, project = "Balboa vivo 1 month", min.cells = 3, min.features = 200)
original_dataset = rep("Balboa_SC_vivo_1_month_sample_2", 994)
balboa_vivo_1_2$original_dataset = original_dataset

balboa_vivo_1_data_3 = Read10X(data.dir = "")
balboa_vivo_1_3 = CreateSeuratObject(balboa_vivo_1_data_3, project = "Balboa vivo 1 month", min.cells = 3, min.features = 200)
original_dataset = rep("Balboa_SC_vivo_1_month_sample_3", 1619)
balboa_vivo_1_3$original_dataset = original_dataset

balboa_vivo_1 = merge(balboa_vivo_1_1, c(balboa_vivo_1_2,balboa_vivo_1_3 ))

## Quality Control & Removal of Mouse Cells
VlnPlot(balboa_vivo_1, features = "TTC36")
TTC36_selection = GetAssayData(object = balboa_vivo_1, 
                               assay = "RNA", slot = "data")["TTC36",]
first_filter = names(which(TTC36_selection<0.5))
balboa_vivo_1 = subset(balboa_vivo_1,cells=first_filter)

balboa_vivo_1[["percent.mt"]] = PercentageFeatureSet(balboa_vivo_1, pattern = "^MT-")
VlnPlot(balboa_vivo_1, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                                 axis.title.y = element_text(size=0),
                                                                                                 axis.text.x = element_text(size=0),
                                                                                                 axis.text.y = element_text(size=20),
                                                                                                 title = element_text(size=12))
balboa_vivo_1_QC = subset(balboa_vivo_1, subset = percent.mt < 25 & nFeature_RNA > 350 & nCount_RNA < 20000)

## Adding Metadata
age = rep("SC-Vivo", 3477)
BMI = rep("SC-Vivo", 3477)
gender = rep(NA, 3477)
balboa_vivo_1_QC$age = age
balboa_vivo_1_QC$BMI = BMI
balboa_vivo_1_QC$gender = gender

## PCA & Clustering
balboa_vivo_1_QC <- NormalizeData(balboa_vivo_1_QC, normalization.method = "LogNormalize", scale.factor = 10000)
balboa_vivo_1_QC <- FindVariableFeatures(balboa_vivo_1_QC, selection.method = "vst", nfeatures = 2000)
balboa_vivo_1_QC <- ScaleData(balboa_vivo_1_QC)
balboa_vivo_1_QC <- RunPCA(balboa_vivo_1_QC, features = VariableFeatures(object = balboa_vivo_1_QC))
balboa_vivo_1_QC <- FindNeighbors(balboa_vivo_1_QC, dims = 1:20)
balboa_vivo_1_QC <- FindClusters(balboa_vivo_1_QC, resolution = 4.5)
balboa_vivo_1_QC <- RunUMAP(balboa_vivo_1_QC, dims = 1:20)
DimPlot(balboa_vivo_1_QC, reduction = "umap", label = T)
FeaturePlot(balboa_vivo_1_QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
balboa_vivo_1_QC = RenameIdents(balboa_vivo_1_QC, '0' = "Beta", '1' = "Alpha",
                                '2' = "Alpha", '3' = "CHGA-", '4' = "Beta", 
                                '5' = "CHGA-", '6'="Alpha", '7'="Alpha",
                                '8'= "EC",'9'="Alpha", '10'="CHGA-", '11'="CHGA-",
                                '12'= "EC",'13' = "CHGA-", '14'="CHGA-", '15'="Beta",
                                '16'= "CHGA-",'17'= "EC", '18' = "CHGA-", '19' = "CHGA-", '20' = "CHGA-",
                                '21' = "CHGA-", '22' = "CHGA-", '23' = "CHGA-", '24' = "CHGA-",
                                '25' = "CHGA-", '26' = "Delta", '27' = "CHGA-", '28' = "Alpha", 
                                '29' = "CHGA-", '30' = "EC",'31' = "EC", '32' = "Prolif", 
                                '33' = "Beta")
balboa_vivo_1_QC$first_identity = balboa_vivo_1_QC@active.ident

## UMAP
DimPlot(balboa_vivo_1_QC, reduction = "umap", label = F, cols = c("#E41A1C", "#377EB8", "grey", "#4DAF4A", "#984EA3", "#7570B3"))+NoLegend()

## Subsetting CHGA+ Cells
balboa_vivo_1_endocrine = subset(balboa_vivo_1_QC, identity = c('Beta', 'Alpha', 'EC', 'Delta', 'Prolif'))


# Downloading H1 6-month TXP SC-islet raw data from Balboa et al. 2022 (DOI:10.1038/s41587-022-01219-z)
# GSE167880 (GSM: GSM5114476)

## Download Raw Data & Creating Seurat Object
setwd("")
balboa_vivo_6_data = Read10X(data.dir = "")
balboa_vivo_6 = CreateSeuratObject(balboa_vivo_6_data, project = "Balboa vivo 6 month", min.cells = 3, min.features = 200)

## Quality Control & Removal of Mouse Cells
VlnPlot(balboa_vivo_6, features = "TTC36")
TTC36_selection = GetAssayData(object = balboa_vivo_6, 
                               assay = "RNA", slot = "data")["TTC36",]
first_filter = names(which(TTC36_selection<0.5))
balboa_vivo_6 = subset(balboa_vivo_6,cells=first_filter)

balboa_vivo_6[["percent.mt"]] = PercentageFeatureSet(balboa_vivo_6, pattern = "^MT-")
VlnPlot(balboa_vivo_6, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                                 axis.title.y = element_text(size=0),
                                                                                                 axis.text.x = element_text(size=0),
                                                                                                 axis.text.y = element_text(size=20),
                                                                                                 title = element_text(size=12))
balboa_vivo_6_QC = subset(balboa_vivo_6, subset = percent.mt < 25 & nFeature_RNA > 1250 & nCount_RNA < 50000)

## Adding Metadata
original_dataset = rep("Balboa_SC_vivo_6_month_sample_1", 2465)
balboa_vivo_6_QC$original_dataset = original_dataset
age = rep("SC-Vivo", 2465)
BMI = rep("SC-Vivo", 2465)
gender = rep(NA, 2465)
balboa_vivo_6_QC$age = age
balboa_vivo_6_QC$BMI = BMI
balboa_vivo_6_QC$gender = gender

## PCA & Clustering
balboa_vivo_6_QC <- NormalizeData(balboa_vivo_6_QC, normalization.method = "LogNormalize", scale.factor = 10000)
balboa_vivo_6_QC <- FindVariableFeatures(balboa_vivo_6_QC, selection.method = "vst", nfeatures = 2000)
balboa_vivo_6_QC <- ScaleData(balboa_vivo_6_QC)
balboa_vivo_6_QC <- RunPCA(balboa_vivo_6_QC, features = VariableFeatures(object = balboa_vivo_6_QC))
balboa_vivo_6_QC <- FindNeighbors(balboa_vivo_6_QC, dims = 1:20)
balboa_vivo_6_QC <- FindClusters(balboa_vivo_6_QC, resolution = 1)
balboa_vivo_6_QC <- RunUMAP(balboa_vivo_6_QC, dims = 1:20)
DimPlot(balboa_vivo_6_QC, reduction = "umap", label = T)
FeaturePlot(balboa_vivo_6_QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
balboa_vivo_6_QC = RenameIdents(balboa_vivo_6_QC, '0' = "Alpha", '1' = "CHGA-",
                                '2' = "Alpha", '3' = "CHGA-", '4' = "Beta", 
                                '5' = "Alpha", '6'="Alpha", '7'="Alpha",
                                '8'= "CHGA-",'9'="Beta", '10'="EC", '11'="CHGA-",
                                '12'= "Beta",'13' = "EC", '14'="CHGA-")
balboa_vivo_6_QC$first_identity = balboa_vivo_6_QC@active.ident

## UMAP
DimPlot(balboa_vivo_6_QC, reduction = "umap", label = F, cols = c("#377EB8", "grey", "#E41A1C", "#4DAF4A"))+NoLegend()

## Subsetting CHGA+ Cells
balboa_vivo_6_endocrine = subset(balboa_vivo_6_QC, identity = c('Alpha', 'Beta', 'EC'))

# Downloading H1 6-month TXP SC-islet raw data from Augsornworawat et al. 2020 (DOI:10.1016/j.celrep.2020.108067)
# GSE151117 (GSM: GSM4567001, GSM4567002, & GSM4567003)

## Download Raw Data, Editing Gene Names, & Creating Seurat Object
setwd("")
augs_vivo_data_1 = Read10X(data.dir = "")
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5I"] ="ATP5ME"
rownames(augs_vivo_data_1)[colnames(augs_vivo_data_1) == "ATP5G2"] ="ATP5MC2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5E"] ="ATP5F1E"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5L"] ="ATP5MG"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "TCEB2"] ="ELOB"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SHFM1"] ="SEM1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "USMG5"] ="ATP5MD"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "NGFRAP1"] ="BEX3"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM213A"] ="PRXL2A"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "GNB2L1"] ="RACK1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "TCEB3"] ="ELOA"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C8orf59"] ="RBIS"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C19orf24"] ="FAM174C"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C20orf24"] ="RAB5IF"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM173A"] ="ANTKMT"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C6orf1"] ="SMIM29"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SEPT2"] ="SEPTIN2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM57B"] ="TLCD3B"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "PVRL2"] ="NECTIN2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM96A"] ="CIAO2A"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "KIAA1715"] ="LNPK"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "RTFDC1"] ="RTF2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM134A"] ="RETREG2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "STRA13"] ="CENPX"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C9orf3"] ="AOPEP"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "HBT8"] ="PWAR6"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "PVRL3"] ="NECTIN3"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "WRB"] ="GET1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ASNA1"] ="GET3"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM60A"] ="SINHCAF"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM63B"] ="MINDY2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "COL4A3BP"] ="CERT1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "WHSC1"] ="NSD2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "APOPT1"] ="COA8"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "TOMM70A"] ="TOMM70"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "PVRL1"] ="NECTIN1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "TMEM261"] ="DMAC1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM195A"] ="MCRIP2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "PPP2R4"] ="PTPA"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SUMO2.1"] ="SUMO2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "UQR11.1"] ="UQCR11"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C11orf73"] ="HIKESHI"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C19orf60"] ="REX1BD"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ZCCHC11"] ="TUT4"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "LHFP"] ="LHFPL6"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "MESDC2"] ="MESD"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM127B"] ="RTL8A"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C7orf55"] ="FMC1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "LINC00493"] ="SMIM26"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "WBSCR22"] ="BUD23"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM105A"] ="OTULINL"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM127A"] ="RTL8C"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C14orf1"] ="ERG28"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SEPT2"] ="SEPTIN2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "VIMP"] ="SELENOS"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C16orf13"] ="METTL26"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SELM"] ="SELENOM"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM195B"] ="MCRIP1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "MLLT4"] ="AFDN"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "WBP5"] ="TCEAL9"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SEPT11"] ="SEPTIN11"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C19orf70"] ="MICOS13"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "WHSC1L1"] ="NSD3"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5H"] ="ATP5PD"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "TCEB1"] ="ELOC"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5C1"] ="ATP5F1C"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SEP15"] ="SELENOF"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM96B"] ="CIA02B"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C17orf89"] ="NDUFAF8"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "APOA1BP"] ="NAXE"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5F1"] ="ATP5PB"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM46A"] ="TENT5A"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C14orf166"] ="RTRAF"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SEPT7"] ="SEPTIN7"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "MINOS1"] ="MICOS10"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C11orf31"] ="SELENOH"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "GLTSCR2"] ="NOP53"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "LINC01420"] ="NBDY"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "FAM159B"] ="SHISAL2B"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C7orf73"] ="STMP1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C19orf43"] ="TRIR"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SELT"] ="SELENOT"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "AES"] ="TLE5"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5G1"] ="ATP5MC1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5B"] ="ATP5F1B"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SEPW1"] ="SELENOW"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5D"] ="ATP5F1D"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5O"] ="ATP5PO"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SELK"] ="SELENOK"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "HN1"] ="JPT1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATPIF1"] ="ATP5IF1"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5J"] ="ATP5PF"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5A1"] ="ATP5F1A"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5G3"] ="ATP5MC3"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "ATP5J2"] ="ATP5MF"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "MYEOV2"] ="COPS9"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "C14orf2"] ="ATP5MPL"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "UQCR11.1"] ="UQCR11"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SUMO2.1"] ="SUMO2"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "LRRC75A-AS1"] ="SNHG29"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "SEPT4"] ="SEPTIN4"
rownames(augs_vivo_data_1)[rownames(augs_vivo_data_1) == "LARGE1"] ="LARGE"
augs_vivo_1 = CreateSeuratObject(augs_vivo_data_1, project = "Augs vivo", min.cells = 3, min.features = 200)
original_dataset = rep("Aug_SC_vivo_6_month_sample_1", 4003)
augs_vivo_1$original_dataset = original_dataset

augs_vivo_data_2 = Read10X(data.dir = "")
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5I"] ="ATP5ME"
rownames(augs_vivo_data_2)[colnames(augs_vivo_data_2) == "ATP5G2"] ="ATP5MC2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5E"] ="ATP5F1E"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5L"] ="ATP5MG"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "TCEB2"] ="ELOB"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SHFM1"] ="SEM1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "USMG5"] ="ATP5MD"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "NGFRAP1"] ="BEX3"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM213A"] ="PRXL2A"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "GNB2L1"] ="RACK1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "TCEB3"] ="ELOA"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C8orf59"] ="RBIS"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C19orf24"] ="FAM174C"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C20orf24"] ="RAB5IF"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM173A"] ="ANTKMT"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C6orf1"] ="SMIM29"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SEPT2"] ="SEPTIN2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM57B"] ="TLCD3B"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "PVRL2"] ="NECTIN2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM96A"] ="CIAO2A"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "KIAA1715"] ="LNPK"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "RTFDC1"] ="RTF2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM134A"] ="RETREG2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "STRA13"] ="CENPX"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C9orf3"] ="AOPEP"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "HBT8"] ="PWAR6"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "PVRL3"] ="NECTIN3"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "WRB"] ="GET1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ASNA1"] ="GET3"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM60A"] ="SINHCAF"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM63B"] ="MINDY2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "COL4A3BP"] ="CERT1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "WHSC1"] ="NSD2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "APOPT1"] ="COA8"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "TOMM70A"] ="TOMM70"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "PVRL1"] ="NECTIN1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "TMEM261"] ="DMAC1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM195A"] ="MCRIP2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "PPP2R4"] ="PTPA"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SUMO2.1"] ="SUMO2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "UQR11.1"] ="UQCR11"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C11orf73"] ="HIKESHI"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C19orf60"] ="REX1BD"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ZCCHC11"] ="TUT4"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "LHFP"] ="LHFPL6"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "MESDC2"] ="MESD"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM127B"] ="RTL8A"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C7orf55"] ="FMC1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "LINC00493"] ="SMIM26"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "WBSCR22"] ="BUD23"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM105A"] ="OTULINL"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM127A"] ="RTL8C"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C14orf1"] ="ERG28"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SEPT2"] ="SEPTIN2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "VIMP"] ="SELENOS"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C16orf13"] ="METTL26"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SELM"] ="SELENOM"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM195B"] ="MCRIP1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "MLLT4"] ="AFDN"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "WBP5"] ="TCEAL9"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SEPT11"] ="SEPTIN11"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C19orf70"] ="MICOS13"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "WHSC1L1"] ="NSD3"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5H"] ="ATP5PD"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "TCEB1"] ="ELOC"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5C1"] ="ATP5F1C"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SEP15"] ="SELENOF"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM96B"] ="CIA02B"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C17orf89"] ="NDUFAF8"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "APOA1BP"] ="NAXE"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5F1"] ="ATP5PB"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM46A"] ="TENT5A"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C14orf166"] ="RTRAF"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SEPT7"] ="SEPTIN7"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "MINOS1"] ="MICOS10"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C11orf31"] ="SELENOH"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "GLTSCR2"] ="NOP53"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "LINC01420"] ="NBDY"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "FAM159B"] ="SHISAL2B"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C7orf73"] ="STMP1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C19orf43"] ="TRIR"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SELT"] ="SELENOT"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "AES"] ="TLE5"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5G1"] ="ATP5MC1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5B"] ="ATP5F1B"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SEPW1"] ="SELENOW"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5D"] ="ATP5F1D"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5O"] ="ATP5PO"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SELK"] ="SELENOK"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "HN1"] ="JPT1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATPIF1"] ="ATP5IF1"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5J"] ="ATP5PF"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5A1"] ="ATP5F1A"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5G3"] ="ATP5MC3"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "ATP5J2"] ="ATP5MF"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "MYEOV2"] ="COPS9"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "C14orf2"] ="ATP5MPL"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "UQCR11.1"] ="UQCR11"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SUMO2.1"] ="SUMO2"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "LRRC75A-AS1"] ="SNHG29"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "SEPT4"] ="SEPTIN4"
rownames(augs_vivo_data_2)[rownames(augs_vivo_data_2) == "LARGE1"] ="LARGE"
augs_vivo_2 = CreateSeuratObject(augs_vivo_data_2, project = "Augs vivo", min.cells = 3, min.features = 200)
original_dataset = rep("Aug_SC_vivo_6_month_sample_2", 2616)
augs_vivo_2$original_dataset = original_dataset

augs_vivo_data_3 = Read10X(data.dir = "")
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5I"] ="ATP5ME"
rownames(augs_vivo_data_3)[colnames(augs_vivo_data_3) == "ATP5G2"] ="ATP5MC2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5E"] ="ATP5F1E"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5L"] ="ATP5MG"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "TCEB2"] ="ELOB"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SHFM1"] ="SEM1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "USMG5"] ="ATP5MD"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "NGFRAP1"] ="BEX3"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM213A"] ="PRXL2A"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "GNB2L1"] ="RACK1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "TCEB3"] ="ELOA"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C8orf59"] ="RBIS"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C19orf24"] ="FAM174C"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C20orf24"] ="RAB5IF"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM173A"] ="ANTKMT"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C6orf1"] ="SMIM29"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SEPT2"] ="SEPTIN2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM57B"] ="TLCD3B"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "PVRL2"] ="NECTIN2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM96A"] ="CIAO2A"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "KIAA1715"] ="LNPK"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "RTFDC1"] ="RTF2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM134A"] ="RETREG2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "STRA13"] ="CENPX"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C9orf3"] ="AOPEP"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "HBT8"] ="PWAR6"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "PVRL3"] ="NECTIN3"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "WRB"] ="GET1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ASNA1"] ="GET3"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM60A"] ="SINHCAF"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM63B"] ="MINDY2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "COL4A3BP"] ="CERT1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "WHSC1"] ="NSD2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "APOPT1"] ="COA8"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "TOMM70A"] ="TOMM70"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "PVRL1"] ="NECTIN1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "TMEM261"] ="DMAC1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM195A"] ="MCRIP2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "PPP2R4"] ="PTPA"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SUMO2.1"] ="SUMO2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "UQR11.1"] ="UQCR11"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C11orf73"] ="HIKESHI"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C19orf60"] ="REX1BD"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ZCCHC11"] ="TUT4"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "LHFP"] ="LHFPL6"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "MESDC2"] ="MESD"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM127B"] ="RTL8A"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C7orf55"] ="FMC1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "LINC00493"] ="SMIM26"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "WBSCR22"] ="BUD23"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM105A"] ="OTULINL"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM127A"] ="RTL8C"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C14orf1"] ="ERG28"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SEPT2"] ="SEPTIN2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "VIMP"] ="SELENOS"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C16orf13"] ="METTL26"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SELM"] ="SELENOM"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM195B"] ="MCRIP1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "MLLT4"] ="AFDN"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "WBP5"] ="TCEAL9"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SEPT11"] ="SEPTIN11"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C19orf70"] ="MICOS13"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "WHSC1L1"] ="NSD3"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5H"] ="ATP5PD"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "TCEB1"] ="ELOC"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5C1"] ="ATP5F1C"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SEP15"] ="SELENOF"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM96B"] ="CIA02B"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C17orf89"] ="NDUFAF8"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "APOA1BP"] ="NAXE"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5F1"] ="ATP5PB"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM46A"] ="TENT5A"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C14orf166"] ="RTRAF"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SEPT7"] ="SEPTIN7"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "MINOS1"] ="MICOS10"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C11orf31"] ="SELENOH"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "GLTSCR2"] ="NOP53"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "LINC01420"] ="NBDY"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "FAM159B"] ="SHISAL2B"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C7orf73"] ="STMP1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C19orf43"] ="TRIR"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SELT"] ="SELENOT"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "AES"] ="TLE5"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5G1"] ="ATP5MC1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5B"] ="ATP5F1B"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SEPW1"] ="SELENOW"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5D"] ="ATP5F1D"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5O"] ="ATP5PO"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SELK"] ="SELENOK"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "HN1"] ="JPT1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATPIF1"] ="ATP5IF1"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5J"] ="ATP5PF"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5A1"] ="ATP5F1A"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5G3"] ="ATP5MC3"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "ATP5J2"] ="ATP5MF"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "MYEOV2"] ="COPS9"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "C14orf2"] ="ATP5MPL"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "UQCR11.1"] ="UQCR11"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SUMO2.1"] ="SUMO2"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "LRRC75A-AS1"] ="SNHG29"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "SEPT4"] ="SEPTIN4"
rownames(augs_vivo_data_3)[rownames(augs_vivo_data_3) == "LARGE1"] ="LARGE"
augs_vivo_3 = CreateSeuratObject(augs_vivo_data_3, project = "Augs vivo", min.cells = 3, min.features = 200)
original_dataset = rep("Aug_SC_vivo_6_month_sample_3", 9007)
augs_vivo_3$original_dataset = original_dataset

augs_vivo = merge(augs_vivo_1, c(augs_vivo_2, augs_vivo_3))

## Quality Control & Removal of Mouse Cells
TTC36_selection = GetAssayData(object = augs_vivo, 
                               assay = "RNA", slot = "data")["TTC36",]
first_filter = names(which(TTC36_selection<3))
augs_vivo = subset(augs_vivo,cells=first_filter)

augs_vivo[["percent.mt"]] = PercentageFeatureSet(augs_vivo, pattern = "^MT-")
VlnPlot(augs_vivo, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                             axis.title.y = element_text(size=0),
                                                                                             axis.text.x = element_text(size=0),
                                                                                             axis.text.y = element_text(size=20),
                                                                                             title = element_text(size=12))
augs_vivo_QC = subset(augs_vivo, subset = percent.mt < 25 & nFeature_RNA > 400 & nCount_RNA < 20000)

## Adding Metadata
age = rep("SC-Vivo", 2728)
BMI = rep("SC-Vivo", 2728)
gender = rep(NA, 2728)
augs_vivo_QC$age = age
augs_vivo_QC$BMI = BMI
augs_vivo_QC$gender = gender

## PCA & Clustering
augs_vivo_QC <- NormalizeData(augs_vivo_QC, normalization.method = "LogNormalize", scale.factor = 10000)
augs_vivo_QC <- FindVariableFeatures(augs_vivo_QC, selection.method = "vst", nfeatures = 2000)
augs_vivo_QC <- ScaleData(augs_vivo_QC)
augs_vivo_QC <- RunPCA(augs_vivo_QC, features = VariableFeatures(object = augs_vivo_QC))
augs_vivo_QC <- FindNeighbors(augs_vivo_QC, dims = 1:20)
augs_vivo_QC <- FindClusters(augs_vivo_QC, resolution = 0.4)
augs_vivo_QC <- RunUMAP(augs_vivo_QC, dims = 1:20)
DimPlot(augs_vivo_QC, reduction = "umap", label = T)
FeaturePlot(augs_vivo_QC, features = "CHGA", cols = brewer.pal(9, 'OrRd')) 

## Assigning Cell Identity
augs_vivo_QC = RenameIdents(augs_vivo_QC, '0' = "EC", '1' = "Beta",
                            '2' = "CHGA-", '3' = "Alpha", '4' = "CHGA-", 
                            '5' = "CHGA-", '6'="CHGA-", '7'="Beta",
                            '8'= "Exo (CHGA+)",'9'="EC")
augs_vivo_QC$first_identity = augs_vivo_QC@active.ident

## UMAP
DimPlot(augs_vivo_QC, reduction = "umap", label = F, cols = c("#4DAF4A", "#E41A1C", "grey", "#377EB8", "yellow2"))+NoLegend()

## Subsetting CHGA+ Cells
augs_vivo_endocrine = subset(augs_vivo_QC, identity = c('EC', "Beta", "Alpha", 'Exo (CHGA+)'))


# Downloading Adult human islet raw data from Baron et al. 2016 (DOI:10.1016/j.cels.2016.08.011)
# GSE84133 (GSM: GSM2230757, GSM2230758, & GSM2230759)

## Download Raw Data, Editing Gene Names, Creating Seurat Object, & Assigning Metadata
setwd("")
Baron1.counts <- read.csv("")
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX2.1"] ="NKX2-1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX2.2"] ="NKX2-2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX2.3"] ="NKX2-3"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX2.5"] ="NKX2-5"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX2.6"] ="NKX2-6"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX2.8"] ="NKX2-8"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX3.1"] ="NKX3-1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX3.2"] ="NKX3-2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX6.1"] ="NKX6-1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX6.2"] ="NKX6-2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NKX6.3"] ="NKX6-3"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.A"] ="HLA-A"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.B"] ="HLA-B"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.C"] ="HLA-C"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DMA"] ="HLA-DMA"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DMB"] ="HLA-DMB"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DOA"] ="HLA-DOA"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DOB"] ="HLA-DOB"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DPA1"] ="HLA-DPA1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DPB1"] ="HLA-DPB1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DQA1"] ="HLA-DQA1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DQA2"] ="HLA-DQA2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DQB1"] ="HLA-DQB1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DQB2"] ="HLA-DQB2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DRA"] ="HLA-DRA"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DRB1"] ="HLA-DRB1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.DRB5"] ="HLA-DRB5"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.E"] ="HLA-E"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.F"] ="HLA-F"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HLA.G"] ="HLA-G"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5I"] ="ATP5ME"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5G2"] ="ATP5MC2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5E"] ="ATP5F1E"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5L"] ="ATP5MG"
colnames(Baron1.counts)[colnames(Baron1.counts) == "TCEB2"] ="ELOB"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SHFM1"] ="SEM1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "USMG5"] ="ATP5MD"
colnames(Baron1.counts)[colnames(Baron1.counts) == "NGFRAP1"] ="BEX3"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM213A"] ="PRXL2A"
colnames(Baron1.counts)[colnames(Baron1.counts) == "GNB2L1"] ="RACK1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ZMYM6NB"] ="TMEM35B"
colnames(Baron1.counts)[colnames(Baron1.counts) == "TCEB3"] ="ELOA"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C8orf59"] ="RBIS"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C19orf24"] ="FAM174C"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C20orf24"] ="RAB5IF"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM173A"] ="ANTKMT"
colnames(Baron1.counts)[colnames(Baron1.counts) == "PALM2.AKAP2"] ="PALM-AKAP2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C6orf1"] ="SMIM29"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SEPT2"] ="SEPTIN2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM57B"] ="TLCD3B"
colnames(Baron1.counts)[colnames(Baron1.counts) == "PVRL2"] ="NECTIN2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM96A"] ="CIAO2A"
colnames(Baron1.counts)[colnames(Baron1.counts) == "KIAA1715"] ="LNPK"
colnames(Baron1.counts)[colnames(Baron1.counts) == "RTFDC1"] ="RTF2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM134A"] ="RETREG2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "STRA13"] ="CENPX"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C9orf3"] ="AOPEP"
colnames(Baron1.counts)[colnames(Baron1.counts) == "PVRL3"] ="NECTIN3"
colnames(Baron1.counts)[colnames(Baron1.counts) == "WRB"] ="GET1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ASNA1"] ="GET3"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM60A"] ="SINHCAF"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM63B"] ="MINDY2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "COL4A3BP"] ="CERT1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "WHSC1"] ="NSD2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "APOPT1"] ="COA8"
colnames(Baron1.counts)[colnames(Baron1.counts) == "TOMM70A"] ="TOMM70"
colnames(Baron1.counts)[colnames(Baron1.counts) == "PVRL1"] ="NECTIN1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "TMEM261"] ="DMAC1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM195A"] ="MCRIP2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "PPP2R4"] ="PTPA"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C11orf73"] ="HIKESHI"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C19orf60"] ="REX1BD"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ZCCHC11"] ="TUT4"
colnames(Baron1.counts)[colnames(Baron1.counts) == "LHFP"] ="LHFPL6"
colnames(Baron1.counts)[colnames(Baron1.counts) == "MESDC2"] ="MESD"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM127B"] ="RTL8A"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C7orf55"] ="FMC1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "LINC00493"] ="SMIM26"
colnames(Baron1.counts)[colnames(Baron1.counts) == "WBSCR22"] ="BUD23"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM105A"] ="OTULINL"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM127A"] ="RTL8C"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C14orf1"] ="ERG28"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SEPT2"] ="SEPTIN2"
colnames(Baron1.counts)[colnames(Baron1.counts) == "VIMP"] ="SELENOS"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C16orf13"] ="METTL26"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SELM"] ="SELENOM"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM195B"] ="MCRIP1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "MLLT4"] ="AFDN"
colnames(Baron1.counts)[colnames(Baron1.counts) == "WBP5"] ="TCEAL9"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SEPT11"] ="SEPTIN11"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C19orf70"] ="MICOS13"
colnames(Baron1.counts)[colnames(Baron1.counts) == "WHSC1L1"] ="NSD3"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5H"] ="ATP5PD"
colnames(Baron1.counts)[colnames(Baron1.counts) == "TCEB1"] ="ELOC"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5C1"] ="ATP5F1C"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SEP15"] ="SELENOF"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM96B"] ="CIA02B"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C17orf89"] ="NDUFAF8"
colnames(Baron1.counts)[colnames(Baron1.counts) == "APOA1BP"] ="NAXE"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5F1"] ="ATP5PB"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM46A"] ="TENT5A"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C14orf166"] ="RTRAF"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SEPT7"] ="SEPTIN7"
colnames(Baron1.counts)[colnames(Baron1.counts) == "MINOS1"] ="MICOS10"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C11orf31"] ="SELENOH"
colnames(Baron1.counts)[colnames(Baron1.counts) == "GLTSCR2"] ="NOP53"
colnames(Baron1.counts)[colnames(Baron1.counts) == "LINC01420"] ="NBDY"
colnames(Baron1.counts)[colnames(Baron1.counts) == "FAM159B"] ="SHISAL2B"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C7orf73"] ="STMP1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C19orf43"] ="TRIR"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SELT"] ="SELENOT"
colnames(Baron1.counts)[colnames(Baron1.counts) == "AES"] ="TLE5"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5G1"] ="ATP5MC1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5B"] ="ATP5F1B"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SEPW1"] ="SELENOW"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5D"] ="ATP5F1D"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5O"] ="ATP5PO"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SELK"] ="SELENOK"
colnames(Baron1.counts)[colnames(Baron1.counts) == "HN1"] ="JPT1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATPIF1"] ="ATP5IF1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5J"] ="ATP5PF"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5A1"] ="ATP5F1A"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5G3"] ="ATP5MC3"
colnames(Baron1.counts)[colnames(Baron1.counts) == "ATP5J2"] ="ATP5MF"
colnames(Baron1.counts)[colnames(Baron1.counts) == "MYEOV2"] ="COPS9"
colnames(Baron1.counts)[colnames(Baron1.counts) == "C14orf2"] ="ATP5MPL"
colnames(Baron1.counts)[colnames(Baron1.counts) == "MIR7.3HG"] ="MIR7-3HG"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SENP3.EIF4A1"] ="SENP3-EIF4A1"
colnames(Baron1.counts)[colnames(Baron1.counts) == "SEPT4"] ="SEPTIN4"
colnames(Baron1.counts)[colnames(Baron1.counts) == "LARGE1"] ="LARGE"
row.names(Baron1.counts) = make.names(Baron1.counts$barcode, unique=TRUE)
Baron1.meta <- select(Baron1.counts, 3)
names(Baron1.meta)[1] <- "CellType"
Baron1.counts <- select(Baron1.counts, -1,-2,-3)
Baron1.data <- CreateSeuratObject(counts = t(Baron1.counts), project = "Adult-B1", meta.data = Baron1.meta)
Baron1.data$CellType = NULL
original_dataset = rep("Adult_cadaveric_baron_sample_1", 1937)
age = rep(17,1937)
BMI = rep(21.5,1937)
gender = rep("male", 1937)
Baron1.data$original_dataset = original_dataset
Baron1.data$age = age
Baron1.data$BMI = BMI
Baron1.data$gender = gender

Baron2.counts <- read.csv("")
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX2.1"] ="NKX2-1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX2.2"] ="NKX2-2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX2.3"] ="NKX2-3"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX2.5"] ="NKX2-5"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX2.6"] ="NKX2-6"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX2.8"] ="NKX2-8"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX3.1"] ="NKX3-1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX3.2"] ="NKX3-2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX6.1"] ="NKX6-1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX6.2"] ="NKX6-2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NKX6.3"] ="NKX6-3"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.A"] ="HLA-A"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.B"] ="HLA-B"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.C"] ="HLA-C"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DMA"] ="HLA-DMA"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DMB"] ="HLA-DMB"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DOA"] ="HLA-DOA"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DOB"] ="HLA-DOB"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DPA1"] ="HLA-DPA1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DPB1"] ="HLA-DPB1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DQA1"] ="HLA-DQA1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DQA2"] ="HLA-DQA2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DQB1"] ="HLA-DQB1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DQB2"] ="HLA-DQB2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DRA"] ="HLA-DRA"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DRB1"] ="HLA-DRB1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.DRB5"] ="HLA-DRB5"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.E"] ="HLA-E"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.F"] ="HLA-F"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HLA.G"] ="HLA-G"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5I"] ="ATP5ME"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5G2"] ="ATP5MC2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5E"] ="ATP5F1E"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5L"] ="ATP5MG"
colnames(Baron2.counts)[colnames(Baron2.counts) == "TCEB2"] ="ELOB"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SHFM1"] ="SEM1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "USMG5"] ="ATP5MD"
colnames(Baron2.counts)[colnames(Baron2.counts) == "NGFRAP1"] ="BEX3"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM213A"] ="PRXL2A"
colnames(Baron2.counts)[colnames(Baron2.counts) == "GNB2L1"] ="RACK1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ZMYM6NB"] ="TMEM35B"
colnames(Baron2.counts)[colnames(Baron2.counts) == "TCEB3"] ="ELOA"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C8orf59"] ="RBIS"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C19orf24"] ="FAM174C"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C20orf24"] ="RAB5IF"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM173A"] ="ANTKMT"
colnames(Baron2.counts)[colnames(Baron2.counts) == "PALM2.AKAP2"] ="PALM-AKAP2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C6orf1"] ="SMIM29"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SEPT2"] ="SEPTIN2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM57B"] ="TLCD3B"
colnames(Baron2.counts)[colnames(Baron2.counts) == "PVRL2"] ="NECTIN2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM96A"] ="CIAO2A"
colnames(Baron2.counts)[colnames(Baron2.counts) == "KIAA1715"] ="LNPK"
colnames(Baron2.counts)[colnames(Baron2.counts) == "RTFDC1"] ="RTF2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM134A"] ="RETREG2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "STRA13"] ="CENPX"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C9orf3"] ="AOPEP"
colnames(Baron2.counts)[colnames(Baron2.counts) == "PVRL3"] ="NECTIN3"
colnames(Baron2.counts)[colnames(Baron2.counts) == "WRB"] ="GET1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ASNA1"] ="GET3"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM60A"] ="SINHCAF"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM63B"] ="MINDY2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "COL4A3BP"] ="CERT1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "WHSC1"] ="NSD2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "APOPT1"] ="COA8"
colnames(Baron2.counts)[colnames(Baron2.counts) == "TOMM70A"] ="TOMM70"
colnames(Baron2.counts)[colnames(Baron2.counts) == "PVRL1"] ="NECTIN1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "TMEM261"] ="DMAC1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM195A"] ="MCRIP2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "PPP2R4"] ="PTPA"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C11orf73"] ="HIKESHI"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C19orf60"] ="REX1BD"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ZCCHC11"] ="TUT4"
colnames(Baron2.counts)[colnames(Baron2.counts) == "LHFP"] ="LHFPL6"
colnames(Baron2.counts)[colnames(Baron2.counts) == "MESDC2"] ="MESD"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM127B"] ="RTL8A"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C7orf55"] ="FMC1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "LINC00493"] ="SMIM26"
colnames(Baron2.counts)[colnames(Baron2.counts) == "WBSCR22"] ="BUD23"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM105A"] ="OTULINL"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM127A"] ="RTL8C"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C14orf1"] ="ERG28"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SEPT2"] ="SEPTIN2"
colnames(Baron2.counts)[colnames(Baron2.counts) == "VIMP"] ="SELENOS"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C16orf13"] ="METTL26"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SELM"] ="SELENOM"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM195B"] ="MCRIP1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "MLLT4"] ="AFDN"
colnames(Baron2.counts)[colnames(Baron2.counts) == "WBP5"] ="TCEAL9"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SEPT11"] ="SEPTIN11"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C19orf70"] ="MICOS13"
colnames(Baron2.counts)[colnames(Baron2.counts) == "WHSC1L1"] ="NSD3"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5H"] ="ATP5PD"
colnames(Baron2.counts)[colnames(Baron2.counts) == "TCEB1"] ="ELOC"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5C1"] ="ATP5F1C"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SEP15"] ="SELENOF"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM96B"] ="CIA02B"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C17orf89"] ="NDUFAF8"
colnames(Baron2.counts)[colnames(Baron2.counts) == "APOA1BP"] ="NAXE"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5F1"] ="ATP5PB"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM46A"] ="TENT5A"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C14orf166"] ="RTRAF"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SEPT7"] ="SEPTIN7"
colnames(Baron2.counts)[colnames(Baron2.counts) == "MINOS1"] ="MICOS10"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C11orf31"] ="SELENOH"
colnames(Baron2.counts)[colnames(Baron2.counts) == "GLTSCR2"] ="NOP53"
colnames(Baron2.counts)[colnames(Baron2.counts) == "LINC01420"] ="NBDY"
colnames(Baron2.counts)[colnames(Baron2.counts) == "FAM159B"] ="SHISAL2B"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C7orf73"] ="STMP1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C19orf43"] ="TRIR"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SELT"] ="SELENOT"
colnames(Baron2.counts)[colnames(Baron2.counts) == "AES"] ="TLE5"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5G1"] ="ATP5MC1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5B"] ="ATP5F1B"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SEPW1"] ="SELENOW"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5D"] ="ATP5F1D"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5O"] ="ATP5PO"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SELK"] ="SELENOK"
colnames(Baron2.counts)[colnames(Baron2.counts) == "HN1"] ="JPT1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATPIF1"] ="ATP5IF1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5J"] ="ATP5PF"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5A1"] ="ATP5F1A"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5G3"] ="ATP5MC3"
colnames(Baron2.counts)[colnames(Baron2.counts) == "ATP5J2"] ="ATP5MF"
colnames(Baron2.counts)[colnames(Baron2.counts) == "MYEOV2"] ="COPS9"
colnames(Baron2.counts)[colnames(Baron2.counts) == "C14orf2"] ="ATP5MPL"
colnames(Baron2.counts)[colnames(Baron2.counts) == "MIR7.3HG"] ="MIR7-3HG"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SENP3.EIF4A1"] ="SENP3-EIF4A1"
colnames(Baron2.counts)[colnames(Baron2.counts) == "SEPT4"] ="SEPTIN4"
colnames(Baron2.counts)[colnames(Baron2.counts) == "LARGE1"] ="LARGE"
row.names(Baron2.counts) = make.names(Baron2.counts$barcode, unique=TRUE)
Baron2.meta <- select(Baron2.counts, 3)
names(Baron2.meta)[1] <- "CellType"
Baron2.counts <- select(Baron2.counts, -1,-2,-3)
Baron2.data <- CreateSeuratObject(counts = t(Baron2.counts), project = "Adult-B2", meta.data = Baron2.meta)
Baron2.data$CellType = NULL
original_dataset = rep("Adult_cadaveric_baron_sample_2", 1724)
age = rep(51,1724)
BMI = rep(21.1,1724)
gender = rep("female", 1724)
Baron2.data$original_dataset = original_dataset
Baron2.data$age = age
Baron2.data$BMI = BMI
Baron2.data$gender = gender

Baron3.counts <- read.csv("")
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX2.1"] ="NKX2-1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX2.2"] ="NKX2-2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX2.3"] ="NKX2-3"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX2.5"] ="NKX2-5"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX2.6"] ="NKX2-6"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX2.8"] ="NKX2-8"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX3.1"] ="NKX3-1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX3.2"] ="NKX3-2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX6.1"] ="NKX6-1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX6.2"] ="NKX6-2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NKX6.3"] ="NKX6-3"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.A"] ="HLA-A"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.B"] ="HLA-B"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.C"] ="HLA-C"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DMA"] ="HLA-DMA"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DMB"] ="HLA-DMB"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DOA"] ="HLA-DOA"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DOB"] ="HLA-DOB"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DPA1"] ="HLA-DPA1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DPB1"] ="HLA-DPB1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DQA1"] ="HLA-DQA1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DQA2"] ="HLA-DQA2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DQB1"] ="HLA-DQB1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DQB2"] ="HLA-DQB2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DRA"] ="HLA-DRA"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DRB1"] ="HLA-DRB1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.DRB5"] ="HLA-DRB5"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.E"] ="HLA-E"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.F"] ="HLA-F"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HLA.G"] ="HLA-G"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5I"] ="ATP5ME"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5G2"] ="ATP5MC2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5E"] ="ATP5F1E"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5L"] ="ATP5MG"
colnames(Baron3.counts)[colnames(Baron3.counts) == "TCEB2"] ="ELOB"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SHFM1"] ="SEM1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "USMG5"] ="ATP5MD"
colnames(Baron3.counts)[colnames(Baron3.counts) == "NGFRAP1"] ="BEX3"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM213A"] ="PRXL2A"
colnames(Baron3.counts)[colnames(Baron3.counts) == "GNB2L1"] ="RACK1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ZMYM6NB"] ="TMEM35B"
colnames(Baron3.counts)[colnames(Baron3.counts) == "TCEB3"] ="ELOA"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C8orf59"] ="RBIS"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C19orf24"] ="FAM174C"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C20orf24"] ="RAB5IF"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM173A"] ="ANTKMT"
colnames(Baron3.counts)[colnames(Baron3.counts) == "PALM2.AKAP2"] ="PALM-AKAP2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C6orf1"] ="SMIM29"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SEPT2"] ="SEPTIN2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM57B"] ="TLCD3B"
colnames(Baron3.counts)[colnames(Baron3.counts) == "PVRL2"] ="NECTIN2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM96A"] ="CIAO2A"
colnames(Baron3.counts)[colnames(Baron3.counts) == "KIAA1715"] ="LNPK"
colnames(Baron3.counts)[colnames(Baron3.counts) == "RTFDC1"] ="RTF2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM134A"] ="RETREG2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "STRA13"] ="CENPX"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C9orf3"] ="AOPEP"
colnames(Baron3.counts)[colnames(Baron3.counts) == "PVRL3"] ="NECTIN3"
colnames(Baron3.counts)[colnames(Baron3.counts) == "WRB"] ="GET1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ASNA1"] ="GET3"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM60A"] ="SINHCAF"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM63B"] ="MINDY2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "COL4A3BP"] ="CERT1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "WHSC1"] ="NSD2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "APOPT1"] ="COA8"
colnames(Baron3.counts)[colnames(Baron3.counts) == "TOMM70A"] ="TOMM70"
colnames(Baron3.counts)[colnames(Baron3.counts) == "PVRL1"] ="NECTIN1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "TMEM261"] ="DMAC1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM195A"] ="MCRIP2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "PPP2R4"] ="PTPA"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C11orf73"] ="HIKESHI"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C19orf60"] ="REX1BD"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ZCCHC11"] ="TUT4"
colnames(Baron3.counts)[colnames(Baron3.counts) == "LHFP"] ="LHFPL6"
colnames(Baron3.counts)[colnames(Baron3.counts) == "MESDC2"] ="MESD"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM127B"] ="RTL8A"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C7orf55"] ="FMC1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "LINC00493"] ="SMIM26"
colnames(Baron3.counts)[colnames(Baron3.counts) == "WBSCR22"] ="BUD23"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM105A"] ="OTULINL"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM127A"] ="RTL8C"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C14orf1"] ="ERG28"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SEPT2"] ="SEPTIN2"
colnames(Baron3.counts)[colnames(Baron3.counts) == "VIMP"] ="SELENOS"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C16orf13"] ="METTL26"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SELM"] ="SELENOM"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM195B"] ="MCRIP1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "MLLT4"] ="AFDN"
colnames(Baron3.counts)[colnames(Baron3.counts) == "WBP5"] ="TCEAL9"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SEPT11"] ="SEPTIN11"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C19orf70"] ="MICOS13"
colnames(Baron3.counts)[colnames(Baron3.counts) == "WHSC1L1"] ="NSD3"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5H"] ="ATP5PD"
colnames(Baron3.counts)[colnames(Baron3.counts) == "TCEB1"] ="ELOC"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5C1"] ="ATP5F1C"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SEP15"] ="SELENOF"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM96B"] ="CIA02B"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C17orf89"] ="NDUFAF8"
colnames(Baron3.counts)[colnames(Baron3.counts) == "APOA1BP"] ="NAXE"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5F1"] ="ATP5PB"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM46A"] ="TENT5A"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C14orf166"] ="RTRAF"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SEPT7"] ="SEPTIN7"
colnames(Baron3.counts)[colnames(Baron3.counts) == "MINOS1"] ="MICOS10"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C11orf31"] ="SELENOH"
colnames(Baron3.counts)[colnames(Baron3.counts) == "GLTSCR2"] ="NOP53"
colnames(Baron3.counts)[colnames(Baron3.counts) == "LINC01420"] ="NBDY"
colnames(Baron3.counts)[colnames(Baron3.counts) == "FAM159B"] ="SHISAL2B"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C7orf73"] ="STMP1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C19orf43"] ="TRIR"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SELT"] ="SELENOT"
colnames(Baron3.counts)[colnames(Baron3.counts) == "AES"] ="TLE5"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5G1"] ="ATP5MC1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5B"] ="ATP5F1B"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SEPW1"] ="SELENOW"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5D"] ="ATP5F1D"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5O"] ="ATP5PO"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SELK"] ="SELENOK"
colnames(Baron3.counts)[colnames(Baron3.counts) == "HN1"] ="JPT1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATPIF1"] ="ATP5IF1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5J"] ="ATP5PF"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5A1"] ="ATP5F1A"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5G3"] ="ATP5MC3"
colnames(Baron3.counts)[colnames(Baron3.counts) == "ATP5J2"] ="ATP5MF"
colnames(Baron3.counts)[colnames(Baron3.counts) == "MYEOV2"] ="COPS9"
colnames(Baron3.counts)[colnames(Baron3.counts) == "C14orf2"] ="ATP5MPL"
colnames(Baron3.counts)[colnames(Baron3.counts) == "MIR7.3HG"] ="MIR7-3HG"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SENP3.EIF4A1"] ="SENP3-EIF4A1"
colnames(Baron3.counts)[colnames(Baron3.counts) == "SEPT4"] ="SEPTIN4"
colnames(Baron3.counts)[colnames(Baron3.counts) == "LARGE1"] ="LARGE"
row.names(Baron3.counts) = make.names(Baron3.counts$barcode, unique=TRUE)
Baron3.meta <- select(Baron3.counts, 3)
names(Baron3.meta)[1] <- "CellType"
Baron3.counts <- select(Baron3.counts, -1,-2,-3)
Baron3.data <- CreateSeuratObject(counts = t(Baron3.counts), project = "Adult-B3", meta.data = Baron3.meta)
Baron3.data$CellType = NULL
original_dataset = rep("Adult_cadaveric_baron_sample_3", 3605)
age = rep(38,3605)
BMI = rep(27.5,3605)
gender = rep("male", 3605)
Baron3.data$original_dataset = original_dataset
Baron3.data$age = age
Baron3.data$BMI = BMI
Baron3.data$gender = gender

Baron.full <- merge(Baron1.data, y = c(Baron2.data,Baron3.data))

## Quality Control 
Baron.full[["percent.mt"]] = PercentageFeatureSet(Baron.full, pattern = "^MT-")
Idents(Baron.full) = "percent.mt"
VlnPlot(Baron.full, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                              axis.title.y = element_text(size=0),
                                                                                              axis.text.x = element_text(size=0),
                                                                                              axis.text.y = element_text(size=20),
                                                                                              title = element_text(size=12))
Baron.QC <- subset(Baron.full, subset = nFeature_RNA < 3000 & nCount_RNA < 10000)

## Adding Metadata
age = rep("SC-Vivo", 2728)
BMI = rep("SC-Vivo", 2728)
gender = rep(NA, 2728)
augs_vivo_QC$age = age
augs_vivo_QC$BMI = BMI
augs_vivo_QC$gender = gender

## PCA & Clustering
Baron.QC <- NormalizeData(Baron.QC)
Baron.QC <- FindVariableFeatures(Baron.QC, selection.method = "vst", nfeatures = 2000)
Baron.QC <- ScaleData(Baron.QC)
Baron.QC <- RunPCA(Baron.QC, features = VariableFeatures(object = Baron.QC))
Baron.QC <- FindNeighbors(Baron.QC, dims = 1:20)
Baron.QC <- FindClusters(Baron.QC, resolution = 0.7)
Baron.QC <- RunUMAP(Baron.QC, dims = 1:20)
DimPlot(Baron.QC, reduction = "umap", label = T)
FeaturePlot(Baron.QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
Baron.QC = RenameIdents(Baron.QC, '0' = "Alpha", '1' = "Beta",
                        '2' = "Alpha", '3' = "Beta", '4' = "Delta", 
                        '5' = "CHGA-", '6'="Beta", '7'="CHGA-",
                        '8'= "CHGA-",'9'="CHGA-", '10'="CHGA-", '11'="PPY",
                        '12'= "Alpha",'13' = "CHGA-", '14'="Beta", '15'="CHGA-",
                        '16'= "CHGA-",'17'= "Alpha", '18' = "CHGA-")
Baron.QC$first_identity = Baron.QC@active.ident

## UMAP
DimPlot(Baron.QC, reduction = "umap", label = F, cols = c("#377EB8", "#E41A1C", "#984EA3", "grey", "#1B9E77"))+NoLegend()

## Subsetting CHGA+ Cells
baron_endocrine = subset(x = Baron.QC, identity = c('Alpha', 'Beta', 'Delta', 'PPY'))


# Downloading Adult human islet raw data from Fang et al. 2019 (DOI:10.1016/j.celrep.2019.02.043)
# GSE101207 (GSM: GSM2700338, GSM2700339, GSM2700340, GSM2700341, & GSM2863188)

## Download Raw Data, Editing Gene Names, Creating Seurat Object, & Assigning Metadata
setwd("")
Fang1.counts <- read.table("",header=T,sep="\t")
rownames(Fang1.counts)<- Fang1.counts[,1]
Fang1.counts <- Fang1.counts[,-1]
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5I"] ="ATP5ME"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5E"] ="ATP5F1E"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5L"] ="ATP5MG"
rownames(Fang1.counts)[rownames(Fang1.counts) == "TCEB2"] ="ELOB"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SHFM1"] ="SEM1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "USMG5"] ="ATP5MD"
rownames(Fang1.counts)[rownames(Fang1.counts) == "NGFRAP1"] ="BEX3"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM213A"] ="PRXL2A"
rownames(Fang1.counts)[rownames(Fang1.counts) == "GNB2L1"] ="RACK1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Fang1.counts)[rownames(Fang1.counts) == "NUPR2"] ="NUPR1L"
rownames(Fang1.counts)[rownames(Fang1.counts) == "TCEB3"] ="ELOA"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C8orf59"] ="RBIS"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C19orf24"] ="FAM174C"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C20orf24"] ="RAB5IF"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM173A"] ="ANTKMT"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C6orf1"] ="SMIM29"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM57B"] ="TLCD3B"
rownames(Fang1.counts)[rownames(Fang1.counts) == "PVRL2"] ="NECTIN2"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM96A"] ="CIAO2A"
rownames(Fang1.counts)[rownames(Fang1.counts) == "KIAA1715"] ="LNPK"
rownames(Fang1.counts)[rownames(Fang1.counts) == "RTFDC1"] ="RTF2"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM134A"] ="RETREG2"
rownames(Fang1.counts)[rownames(Fang1.counts) == "STRA13"] ="CENPX"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C9orf3"] ="AOPEP"
rownames(Fang1.counts)[rownames(Fang1.counts) == "PVRL3"] ="NECTIN3"
rownames(Fang1.counts)[rownames(Fang1.counts) == "WRB"] ="GET1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ASNA1"] ="GET3"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM60A"] ="SINHCAF"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM63B"] ="MINDY2"
rownames(Fang1.counts)[rownames(Fang1.counts) == "COL4A3BP"] ="CERT1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "WHSC1"] ="NSD2"
rownames(Fang1.counts)[rownames(Fang1.counts) == "APOPT1"] ="COA8"
rownames(Fang1.counts)[rownames(Fang1.counts) == "TOMM70A"] ="TOMM70"
rownames(Fang1.counts)[rownames(Fang1.counts) == "PVRL1"] ="NECTIN1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "TMEM261"] ="DMAC1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM195A"] ="MCRIP2"
rownames(Fang1.counts)[rownames(Fang1.counts) == "PPP2R4"] ="PTPA"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C11orf73"] ="HIKESHI"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C19orf60"] ="REX1BD"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ZCCHC11"] ="TUT4"
rownames(Fang1.counts)[rownames(Fang1.counts) == "LHFP"] ="LHFPL6"
rownames(Fang1.counts)[rownames(Fang1.counts) == "MESDC2"] ="MESD"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM127B"] ="RTL8A"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C7orf55"] ="FMC1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "LINC00493"] ="SMIM26"
rownames(Fang1.counts)[rownames(Fang1.counts) == "WBSCR22"] ="BUD23"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM105A"] ="OTULINL"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM127A"] ="RTL8C"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C14orf1"] ="ERG28"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang1.counts)[rownames(Fang1.counts) == "VIMP"] ="SELENOS"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C16orf13"] ="METTL26"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SELM"] ="SELENOM"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM195B"] ="MCRIP1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "MLLT4"] ="AFDN"
rownames(Fang1.counts)[rownames(Fang1.counts) == "WBP5"] ="TCEAL9"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SEPT11"] ="SEPTIN11"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C19orf70"] ="MICOS13"
rownames(Fang1.counts)[rownames(Fang1.counts) == "WHSC1L1"] ="NSD3"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5H"] ="ATP5PD"
rownames(Fang1.counts)[rownames(Fang1.counts) == "TCEB1"] ="ELOC"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SEP15"] ="SELENOF"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM96B"] ="CIA02B"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C17orf89"] ="NDUFAF8"
rownames(Fang1.counts)[rownames(Fang1.counts) == "APOA1BP"] ="NAXE"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5F1"] ="ATP5PB"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM46A"] ="TENT5A"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C14orf166"] ="RTRAF"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SEPT7"] ="SEPTIN7"
rownames(Fang1.counts)[rownames(Fang1.counts) == "MINOS1"] ="MICOS10"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C11orf31"] ="SELENOH"
rownames(Fang1.counts)[rownames(Fang1.counts) == "GLTSCR2"] ="NOP53"
rownames(Fang1.counts)[rownames(Fang1.counts) == "LINC01420"] ="NBDY"
rownames(Fang1.counts)[rownames(Fang1.counts) == "FAM159B"] ="SHISAL2B"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C7orf73"] ="STMP1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C19orf43"] ="TRIR"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SELT"] ="SELENOT"
rownames(Fang1.counts)[rownames(Fang1.counts) == "AES"] ="TLE5"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5B"] ="ATP5F1B"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SEPW1"] ="SELENOW"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5D"] ="ATP5F1D"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5O"] ="ATP5PO"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SELK"] ="SELENOK"
rownames(Fang1.counts)[rownames(Fang1.counts) == "HN1"] ="JPT1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5J"] ="ATP5PF"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Fang1.counts)[rownames(Fang1.counts) == "ATP5J2"] ="ATP5MF"
rownames(Fang1.counts)[rownames(Fang1.counts) == "MYEOV2"] ="COPS9"
rownames(Fang1.counts)[rownames(Fang1.counts) == "C14orf2"] ="ATP5MPL"
rownames(Fang1.counts)[rownames(Fang1.counts) == "LRRC75A-AS1"] ="SNHG29"
rownames(Fang1.counts)[rownames(Fang1.counts) == "SEPT4"] ="SEPTIN4"
rownames(Fang1.counts)[rownames(Fang1.counts) == "LARGE1"] ="LARGE"
Fang1.data <- CreateSeuratObject(counts = Fang1.counts, project = "Adult-F1")
original_dataset = rep("Adult_cadaveric_fang_sample_1", 4000)
age = rep(27,4000)
BMI = rep(20.6,4000)
gender = rep("male", 4000)
Fang1.data$original_dataset = original_dataset
Fang1.data$age = age
Fang1.data$BMI = BMI
Fang1.data$gender = gender

Fang2.counts <- read.table("",header=T,sep="\t")
rownames(Fang2.counts)<- Fang2.counts[,1]
Fang2.counts <- Fang2.counts[,-1]
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5I"] ="ATP5ME"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5E"] ="ATP5F1E"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5L"] ="ATP5MG"
rownames(Fang2.counts)[rownames(Fang2.counts) == "TCEB2"] ="ELOB"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SHFM1"] ="SEM1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "USMG5"] ="ATP5MD"
rownames(Fang2.counts)[rownames(Fang2.counts) == "NGFRAP1"] ="BEX3"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM213A"] ="PRXL2A"
rownames(Fang2.counts)[rownames(Fang2.counts) == "GNB2L1"] ="RACK1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Fang2.counts)[rownames(Fang2.counts) == "NUPR2"] ="NUPR1L"
rownames(Fang2.counts)[rownames(Fang2.counts) == "TCEB3"] ="ELOA"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C8orf59"] ="RBIS"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C19orf24"] ="FAM174C"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C20orf24"] ="RAB5IF"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM173A"] ="ANTKMT"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C6orf1"] ="SMIM29"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM57B"] ="TLCD3B"
rownames(Fang2.counts)[rownames(Fang2.counts) == "PVRL2"] ="NECTIN2"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM96A"] ="CIAO2A"
rownames(Fang2.counts)[rownames(Fang2.counts) == "KIAA1715"] ="LNPK"
rownames(Fang2.counts)[rownames(Fang2.counts) == "RTFDC1"] ="RTF2"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM134A"] ="RETREG2"
rownames(Fang2.counts)[rownames(Fang2.counts) == "STRA13"] ="CENPX"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C9orf3"] ="AOPEP"
rownames(Fang2.counts)[rownames(Fang2.counts) == "PVRL3"] ="NECTIN3"
rownames(Fang2.counts)[rownames(Fang2.counts) == "WRB"] ="GET1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ASNA1"] ="GET3"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM60A"] ="SINHCAF"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM63B"] ="MINDY2"
rownames(Fang2.counts)[rownames(Fang2.counts) == "COL4A3BP"] ="CERT1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "WHSC1"] ="NSD2"
rownames(Fang2.counts)[rownames(Fang2.counts) == "APOPT1"] ="COA8"
rownames(Fang2.counts)[rownames(Fang2.counts) == "TOMM70A"] ="TOMM70"
rownames(Fang2.counts)[rownames(Fang2.counts) == "PVRL1"] ="NECTIN1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "TMEM261"] ="DMAC1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM195A"] ="MCRIP2"
rownames(Fang2.counts)[rownames(Fang2.counts) == "PPP2R4"] ="PTPA"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C11orf73"] ="HIKESHI"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C19orf60"] ="REX1BD"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ZCCHC11"] ="TUT4"
rownames(Fang2.counts)[rownames(Fang2.counts) == "LHFP"] ="LHFPL6"
rownames(Fang2.counts)[rownames(Fang2.counts) == "MESDC2"] ="MESD"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM127B"] ="RTL8A"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C7orf55"] ="FMC1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "LINC00493"] ="SMIM26"
rownames(Fang2.counts)[rownames(Fang2.counts) == "WBSCR22"] ="BUD23"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM105A"] ="OTULINL"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM127A"] ="RTL8C"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C14orf1"] ="ERG28"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang2.counts)[rownames(Fang2.counts) == "VIMP"] ="SELENOS"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C16orf13"] ="METTL26"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SELM"] ="SELENOM"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM195B"] ="MCRIP1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "MLLT4"] ="AFDN"
rownames(Fang2.counts)[rownames(Fang2.counts) == "WBP5"] ="TCEAL9"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SEPT11"] ="SEPTIN11"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C19orf70"] ="MICOS13"
rownames(Fang2.counts)[rownames(Fang2.counts) == "WHSC1L1"] ="NSD3"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5H"] ="ATP5PD"
rownames(Fang2.counts)[rownames(Fang2.counts) == "TCEB1"] ="ELOC"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SEP15"] ="SELENOF"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM96B"] ="CIA02B"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C17orf89"] ="NDUFAF8"
rownames(Fang2.counts)[rownames(Fang2.counts) == "APOA1BP"] ="NAXE"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5F1"] ="ATP5PB"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM46A"] ="TENT5A"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C14orf166"] ="RTRAF"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SEPT7"] ="SEPTIN7"
rownames(Fang2.counts)[rownames(Fang2.counts) == "MINOS1"] ="MICOS10"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C11orf31"] ="SELENOH"
rownames(Fang2.counts)[rownames(Fang2.counts) == "GLTSCR2"] ="NOP53"
rownames(Fang2.counts)[rownames(Fang2.counts) == "LINC01420"] ="NBDY"
rownames(Fang2.counts)[rownames(Fang2.counts) == "FAM159B"] ="SHISAL2B"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C7orf73"] ="STMP1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C19orf43"] ="TRIR"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SELT"] ="SELENOT"
rownames(Fang2.counts)[rownames(Fang2.counts) == "AES"] ="TLE5"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5B"] ="ATP5F1B"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SEPW1"] ="SELENOW"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5D"] ="ATP5F1D"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5O"] ="ATP5PO"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SELK"] ="SELENOK"
rownames(Fang2.counts)[rownames(Fang2.counts) == "HN1"] ="JPT1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5J"] ="ATP5PF"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Fang2.counts)[rownames(Fang2.counts) == "ATP5J2"] ="ATP5MF"
rownames(Fang2.counts)[rownames(Fang2.counts) == "MYEOV2"] ="COPS9"
rownames(Fang2.counts)[rownames(Fang2.counts) == "C14orf2"] ="ATP5MPL"
rownames(Fang2.counts)[rownames(Fang2.counts) == "LRRC75A-AS1"] ="SNHG29"
rownames(Fang2.counts)[rownames(Fang2.counts) == "SEPT4"] ="SEPTIN4"
rownames(Fang2.counts)[rownames(Fang2.counts) == "LARGE1"] ="LARGE"
Fang2.data <- CreateSeuratObject(counts = Fang2.counts, project = "Adult-F2")
original_dataset = rep("Adult_cadaveric_fang_sample_2", 4000)
age = rep(21,4000)
BMI = rep(22.8,4000)
gender = rep("male", 4000)
Fang2.data$original_dataset = original_dataset
Fang2.data$age = age
Fang2.data$BMI = BMI
Fang2.data$gender = gender


Fang3.counts <- read.table("",header=T,sep="\t")
rownames(Fang3.counts)<- Fang3.counts[,1]
Fang3.counts <- Fang3.counts[,-1]
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5I"] ="ATP5ME"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5E"] ="ATP5F1E"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5L"] ="ATP5MG"
rownames(Fang3.counts)[rownames(Fang3.counts) == "TCEB2"] ="ELOB"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SHFM1"] ="SEM1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "USMG5"] ="ATP5MD"
rownames(Fang3.counts)[rownames(Fang3.counts) == "NGFRAP1"] ="BEX3"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM213A"] ="PRXL2A"
rownames(Fang3.counts)[rownames(Fang3.counts) == "GNB2L1"] ="RACK1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Fang3.counts)[rownames(Fang3.counts) == "NUPR2"] ="NUPR1L"
rownames(Fang3.counts)[rownames(Fang3.counts) == "TCEB3"] ="ELOA"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C8orf59"] ="RBIS"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C19orf24"] ="FAM174C"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C20orf24"] ="RAB5IF"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM173A"] ="ANTKMT"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C6orf1"] ="SMIM29"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM57B"] ="TLCD3B"
rownames(Fang3.counts)[rownames(Fang3.counts) == "PVRL2"] ="NECTIN2"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM96A"] ="CIAO2A"
rownames(Fang3.counts)[rownames(Fang3.counts) == "KIAA1715"] ="LNPK"
rownames(Fang3.counts)[rownames(Fang3.counts) == "RTFDC1"] ="RTF2"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM134A"] ="RETREG2"
rownames(Fang3.counts)[rownames(Fang3.counts) == "STRA13"] ="CENPX"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C9orf3"] ="AOPEP"
rownames(Fang3.counts)[rownames(Fang3.counts) == "PVRL3"] ="NECTIN3"
rownames(Fang3.counts)[rownames(Fang3.counts) == "WRB"] ="GET1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ASNA1"] ="GET3"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM60A"] ="SINHCAF"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM63B"] ="MINDY2"
rownames(Fang3.counts)[rownames(Fang3.counts) == "COL4A3BP"] ="CERT1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "WHSC1"] ="NSD2"
rownames(Fang3.counts)[rownames(Fang3.counts) == "APOPT1"] ="COA8"
rownames(Fang3.counts)[rownames(Fang3.counts) == "TOMM70A"] ="TOMM70"
rownames(Fang3.counts)[rownames(Fang3.counts) == "PVRL1"] ="NECTIN1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "TMEM261"] ="DMAC1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM195A"] ="MCRIP2"
rownames(Fang3.counts)[rownames(Fang3.counts) == "PPP2R4"] ="PTPA"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C11orf73"] ="HIKESHI"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C19orf60"] ="REX1BD"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ZCCHC11"] ="TUT4"
rownames(Fang3.counts)[rownames(Fang3.counts) == "LHFP"] ="LHFPL6"
rownames(Fang3.counts)[rownames(Fang3.counts) == "MESDC2"] ="MESD"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM127B"] ="RTL8A"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C7orf55"] ="FMC1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "LINC00493"] ="SMIM26"
rownames(Fang3.counts)[rownames(Fang3.counts) == "WBSCR22"] ="BUD23"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM105A"] ="OTULINL"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM127A"] ="RTL8C"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C14orf1"] ="ERG28"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang3.counts)[rownames(Fang3.counts) == "VIMP"] ="SELENOS"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C16orf13"] ="METTL26"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SELM"] ="SELENOM"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM195B"] ="MCRIP1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "MLLT4"] ="AFDN"
rownames(Fang3.counts)[rownames(Fang3.counts) == "WBP5"] ="TCEAL9"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SEPT11"] ="SEPTIN11"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C19orf70"] ="MICOS13"
rownames(Fang3.counts)[rownames(Fang3.counts) == "WHSC1L1"] ="NSD3"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5H"] ="ATP5PD"
rownames(Fang3.counts)[rownames(Fang3.counts) == "TCEB1"] ="ELOC"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SEP15"] ="SELENOF"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM96B"] ="CIA02B"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C17orf89"] ="NDUFAF8"
rownames(Fang3.counts)[rownames(Fang3.counts) == "APOA1BP"] ="NAXE"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5F1"] ="ATP5PB"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM46A"] ="TENT5A"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C14orf166"] ="RTRAF"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SEPT7"] ="SEPTIN7"
rownames(Fang3.counts)[rownames(Fang3.counts) == "MINOS1"] ="MICOS10"
rownames(Fang3.counts)[rownames(Fang3.counts) == "GLTSCR2"] ="NOP53"
rownames(Fang3.counts)[rownames(Fang3.counts) == "LINC01420"] ="NBDY"
rownames(Fang3.counts)[rownames(Fang3.counts) == "FAM159B"] ="SHISAL2B"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C7orf73"] ="STMP1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C19orf43"] ="TRIR"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SELT"] ="SELENOT"
rownames(Fang3.counts)[rownames(Fang3.counts) == "AES"] ="TLE5"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5B"] ="ATP5F1B"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SEPW1"] ="SELENOW"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5D"] ="ATP5F1D"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5O"] ="ATP5PO"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SELK"] ="SELENOK"
rownames(Fang3.counts)[rownames(Fang3.counts) == "HN1"] ="JPT1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5J"] ="ATP5PF"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Fang3.counts)[rownames(Fang3.counts) == "ATP5J2"] ="ATP5MF"
rownames(Fang3.counts)[rownames(Fang3.counts) == "MYEOV2"] ="COPS9"
rownames(Fang3.counts)[rownames(Fang3.counts) == "C14orf2"] ="ATP5MPL"
rownames(Fang3.counts)[rownames(Fang3.counts) == "LRRC75A-AS1"] ="SNHG29"
rownames(Fang3.counts)[rownames(Fang3.counts) == "SEPT4"] ="SEPTIN4"
rownames(Fang3.counts)[rownames(Fang3.counts) == "LARGE1"] ="LARGE"
Fang3.data <- CreateSeuratObject(counts = Fang3.counts, project = "Adult-F3")
original_dataset = rep("Adult_cadaveric_fang_sample_3", 20000)
age = rep(38,20000)
BMI = rep(34.4,20000)
gender = rep("female", 20000)
Fang3.data$original_dataset = original_dataset
Fang3.data$age = age
Fang3.data$BMI = BMI
Fang3.data$gender = gender

Fang4.counts <- read.table("",header=T,sep="\t")
rownames(Fang4.counts)<- Fang4.counts[,1]
Fang4.counts <- Fang4.counts[,-1]
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5I"] ="ATP5ME"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5E"] ="ATP5F1E"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5L"] ="ATP5MG"
rownames(Fang4.counts)[rownames(Fang4.counts) == "TCEB2"] ="ELOB"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SHFM1"] ="SEM1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "USMG5"] ="ATP5MD"
rownames(Fang4.counts)[rownames(Fang4.counts) == "NGFRAP1"] ="BEX3"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM213A"] ="PRXL2A"
rownames(Fang4.counts)[rownames(Fang4.counts) == "GNB2L1"] ="RACK1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Fang4.counts)[rownames(Fang4.counts) == "NUPR2"] ="NUPR1L"
rownames(Fang4.counts)[rownames(Fang4.counts) == "TCEB3"] ="ELOA"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C8orf59"] ="RBIS"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C19orf24"] ="FAM174C"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C20orf24"] ="RAB5IF"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM173A"] ="ANTKMT"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C6orf1"] ="SMIM29"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM57B"] ="TLCD3B"
rownames(Fang4.counts)[rownames(Fang4.counts) == "PVRL2"] ="NECTIN2"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM96A"] ="CIAO2A"
rownames(Fang4.counts)[rownames(Fang4.counts) == "KIAA1715"] ="LNPK"
rownames(Fang4.counts)[rownames(Fang4.counts) == "RTFDC1"] ="RTF2"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM134A"] ="RETREG2"
rownames(Fang4.counts)[rownames(Fang4.counts) == "STRA13"] ="CENPX"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C9orf3"] ="AOPEP"
rownames(Fang4.counts)[rownames(Fang4.counts) == "PVRL3"] ="NECTIN3"
rownames(Fang4.counts)[rownames(Fang4.counts) == "WRB"] ="GET1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ASNA1"] ="GET3"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM60A"] ="SINHCAF"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM63B"] ="MINDY2"
rownames(Fang4.counts)[rownames(Fang4.counts) == "COL4A3BP"] ="CERT1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "WHSC1"] ="NSD2"
rownames(Fang4.counts)[rownames(Fang4.counts) == "APOPT1"] ="COA8"
rownames(Fang4.counts)[rownames(Fang4.counts) == "TOMM70A"] ="TOMM70"
rownames(Fang4.counts)[rownames(Fang4.counts) == "PVRL1"] ="NECTIN1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "TMEM261"] ="DMAC1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM195A"] ="MCRIP2"
rownames(Fang4.counts)[rownames(Fang4.counts) == "PPP2R4"] ="PTPA"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C11orf73"] ="HIKESHI"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C19orf60"] ="REX1BD"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ZCCHC11"] ="TUT4"
rownames(Fang4.counts)[rownames(Fang4.counts) == "LHFP"] ="LHFPL6"
rownames(Fang4.counts)[rownames(Fang4.counts) == "MESDC2"] ="MESD"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM127B"] ="RTL8A"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C7orf55"] ="FMC1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "LINC00493"] ="SMIM26"
rownames(Fang4.counts)[rownames(Fang4.counts) == "WBSCR22"] ="BUD23"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM105A"] ="OTULINL"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM127A"] ="RTL8C"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C14orf1"] ="ERG28"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang4.counts)[rownames(Fang4.counts) == "VIMP"] ="SELENOS"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C16orf13"] ="METTL26"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SELM"] ="SELENOM"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM195B"] ="MCRIP1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "MLLT4"] ="AFDN"
rownames(Fang4.counts)[rownames(Fang4.counts) == "WBP5"] ="TCEAL9"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SEPT11"] ="SEPTIN11"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C19orf70"] ="MICOS13"
rownames(Fang4.counts)[rownames(Fang4.counts) == "WHSC1L1"] ="NSD3"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5H"] ="ATP5PD"
rownames(Fang4.counts)[rownames(Fang4.counts) == "TCEB1"] ="ELOC"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SEP15"] ="SELENOF"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM96B"] ="CIA02B"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C17orf89"] ="NDUFAF8"
rownames(Fang4.counts)[rownames(Fang4.counts) == "APOA1BP"] ="NAXE"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5F1"] ="ATP5PB"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM46A"] ="TENT5A"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C14orf166"] ="RTRAF"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SEPT7"] ="SEPTIN7"
rownames(Fang4.counts)[rownames(Fang4.counts) == "MINOS1"] ="MICOS10"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C11orf31"] ="SELENOH"
rownames(Fang4.counts)[rownames(Fang4.counts) == "GLTSCR2"] ="NOP53"
rownames(Fang4.counts)[rownames(Fang4.counts) == "LINC01420"] ="NBDY"
rownames(Fang4.counts)[rownames(Fang4.counts) == "FAM159B"] ="SHISAL2B"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C7orf73"] ="STMP1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C19orf43"] ="TRIR"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SELT"] ="SELENOT"
rownames(Fang4.counts)[rownames(Fang4.counts) == "AES"] ="TLE5"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5B"] ="ATP5F1B"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SEPW1"] ="SELENOW"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5D"] ="ATP5F1D"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5O"] ="ATP5PO"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SELK"] ="SELENOK"
rownames(Fang4.counts)[rownames(Fang4.counts) == "HN1"] ="JPT1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5J"] ="ATP5PF"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Fang4.counts)[rownames(Fang4.counts) == "ATP5J2"] ="ATP5MF"
rownames(Fang4.counts)[rownames(Fang4.counts) == "MYEOV2"] ="COPS9"
rownames(Fang4.counts)[rownames(Fang4.counts) == "C14orf2"] ="ATP5MPL"
rownames(Fang4.counts)[rownames(Fang4.counts) == "LRRC75A-AS1"] ="SNHG29"
rownames(Fang4.counts)[rownames(Fang4.counts) == "SEPT4"] ="SEPTIN4"
rownames(Fang4.counts)[rownames(Fang4.counts) == "LARGE1"] ="LARGE"
Fang4.data <- CreateSeuratObject(counts = Fang4.counts, project = "Adult-F4")
original_dataset = rep("Adult_cadaveric_fang_sample_4", 20000)
age = rep(52,20000)
BMI = rep(22,20000)
gender = rep("male", 20000)
Fang4.data$original_dataset = original_dataset
Fang4.data$age = age
Fang4.data$BMI = BMI
Fang4.data$gender = gender

Fang6.counts <- read.table("",header=T,sep="\t")
rownames(Fang6.counts)<- Fang6.counts[,1]
Fang6.counts <- Fang6.counts[,-1]
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5I"] ="ATP5ME"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5E"] ="ATP5F1E"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5L"] ="ATP5MG"
rownames(Fang6.counts)[rownames(Fang6.counts) == "TCEB2"] ="ELOB"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SHFM1"] ="SEM1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "USMG5"] ="ATP5MD"
rownames(Fang6.counts)[rownames(Fang6.counts) == "NGFRAP1"] ="BEX3"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM213A"] ="PRXL2A"
rownames(Fang6.counts)[rownames(Fang6.counts) == "GNB2L1"] ="RACK1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Fang6.counts)[rownames(Fang6.counts) == "NUPR2"] ="NUPR1L"
rownames(Fang6.counts)[rownames(Fang6.counts) == "TCEB3"] ="ELOA"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C8orf59"] ="RBIS"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C19orf24"] ="FAM174C"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C20orf24"] ="RAB5IF"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM173A"] ="ANTKMT"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C6orf1"] ="SMIM29"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM57B"] ="TLCD3B"
rownames(Fang6.counts)[rownames(Fang6.counts) == "PVRL2"] ="NECTIN2"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM96A"] ="CIAO2A"
rownames(Fang6.counts)[rownames(Fang6.counts) == "KIAA1715"] ="LNPK"
rownames(Fang6.counts)[rownames(Fang6.counts) == "RTFDC1"] ="RTF2"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM134A"] ="RETREG2"
rownames(Fang6.counts)[rownames(Fang6.counts) == "STRA13"] ="CENPX"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C9orf3"] ="AOPEP"
rownames(Fang6.counts)[rownames(Fang6.counts) == "PVRL3"] ="NECTIN3"
rownames(Fang6.counts)[rownames(Fang6.counts) == "WRB"] ="GET1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ASNA1"] ="GET3"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM60A"] ="SINHCAF"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM63B"] ="MINDY2"
rownames(Fang6.counts)[rownames(Fang6.counts) == "COL4A3BP"] ="CERT1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "WHSC1"] ="NSD2"
rownames(Fang6.counts)[rownames(Fang6.counts) == "APOPT1"] ="COA8"
rownames(Fang6.counts)[rownames(Fang6.counts) == "TOMM70A"] ="TOMM70"
rownames(Fang6.counts)[rownames(Fang6.counts) == "PVRL1"] ="NECTIN1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "TMEM261"] ="DMAC1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM195A"] ="MCRIP2"
rownames(Fang6.counts)[rownames(Fang6.counts) == "PPP2R4"] ="PTPA"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C11orf73"] ="HIKESHI"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C19orf60"] ="REX1BD"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ZCCHC11"] ="TUT4"
rownames(Fang6.counts)[rownames(Fang6.counts) == "LHFP"] ="LHFPL6"
rownames(Fang6.counts)[rownames(Fang6.counts) == "MESDC2"] ="MESD"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM127B"] ="RTL8A"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C7orf55"] ="FMC1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "LINC00493"] ="SMIM26"
rownames(Fang6.counts)[rownames(Fang6.counts) == "WBSCR22"] ="BUD23"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM105A"] ="OTULINL"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM127A"] ="RTL8C"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C14orf1"] ="ERG28"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SEPT2"] ="SEPTIN2"
rownames(Fang6.counts)[rownames(Fang6.counts) == "VIMP"] ="SELENOS"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C16orf13"] ="METTL26"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SELM"] ="SELENOM"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM195B"] ="MCRIP1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "MLLT4"] ="AFDN"
rownames(Fang6.counts)[rownames(Fang6.counts) == "WBP5"] ="TCEAL9"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SEPT11"] ="SEPTIN11"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C19orf70"] ="MICOS13"
rownames(Fang6.counts)[rownames(Fang6.counts) == "WHSC1L1"] ="NSD3"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5H"] ="ATP5PD"
rownames(Fang6.counts)[rownames(Fang6.counts) == "TCEB1"] ="ELOC"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SEP15"] ="SELENOF"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM96B"] ="CIA02B"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C17orf89"] ="NDUFAF8"
rownames(Fang6.counts)[rownames(Fang6.counts) == "APOA1BP"] ="NAXE"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5F1"] ="ATP5PB"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM46A"] ="TENT5A"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C14orf166"] ="RTRAF"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SEPT7"] ="SEPTIN7"
rownames(Fang6.counts)[rownames(Fang6.counts) == "MINOS1"] ="MICOS10"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C11orf31"] ="SELENOH"
rownames(Fang6.counts)[rownames(Fang6.counts) == "GLTSCR2"] ="NOP53"
rownames(Fang6.counts)[rownames(Fang6.counts) == "LINC01420"] ="NBDY"
rownames(Fang6.counts)[rownames(Fang6.counts) == "FAM159B"] ="SHISAL2B"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C7orf73"] ="STMP1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C19orf43"] ="TRIR"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SELT"] ="SELENOT"
rownames(Fang6.counts)[rownames(Fang6.counts) == "AES"] ="TLE5"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5B"] ="ATP5F1B"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SEPW1"] ="SELENOW"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5D"] ="ATP5F1D"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5O"] ="ATP5PO"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SELK"] ="SELENOK"
rownames(Fang6.counts)[rownames(Fang6.counts) == "HN1"] ="JPT1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5J"] ="ATP5PF"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Fang6.counts)[rownames(Fang6.counts) == "ATP5J2"] ="ATP5MF"
rownames(Fang6.counts)[rownames(Fang6.counts) == "MYEOV2"] ="COPS9"
rownames(Fang6.counts)[rownames(Fang6.counts) == "C14orf2"] ="ATP5MPL"
rownames(Fang6.counts)[rownames(Fang6.counts) == "LRRC75A-AS1"] ="SNHG29"
rownames(Fang6.counts)[rownames(Fang6.counts) == "SEPT4"] ="SEPTIN4"
rownames(Fang6.counts)[rownames(Fang6.counts) == "LARGE1"] ="LARGE"
Fang6.data <- CreateSeuratObject(counts = Fang6.counts, project = "Adult-F6")
original_dataset = rep("Adult_cadaveric_fang_sample_6", 20000)
age = rep(44,20000)
BMI = rep(34.6,20000)
gender = rep("male", 20000)
Fang6.data$original_dataset = original_dataset
Fang6.data$age = age
Fang6.data$BMI = BMI
Fang6.data$gender = gender

Fang.full <- merge(Fang1.data, y = c(Fang2.data,Fang3.data,Fang4.data,Fang6.data), add.cell.ids = c("F1","F2","F3","F4","F6"))


## Quality Control 
Fang.full[["percent.mt"]] = 0
Idents(Fang.full) = "percent.mt"
VlnPlot(Fang.full, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                             axis.title.y = element_text(size=0),
                                                                                             axis.text.x = element_text(size=0),
                                                                                             axis.text.y = element_text(size=20),
                                                                                             title = element_text(size=12))
Fang.QC <- subset(Fang.full, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & nCount_RNA < 5000)

## PCA & Clustering
Fang.QC <- NormalizeData(Fang.QC)
Fang.QC <- FindVariableFeatures(Fang.QC, selection.method = "vst", nfeatures = 2000)
Fang.QC <- ScaleData(Fang.QC)
Fang.QC <- RunPCA(Fang.QC, features = VariableFeatures(object = Fang.QC))
Fang.QC <- FindNeighbors(Fang.QC, dims = 1:20)
Fang.QC <- FindClusters(Fang.QC, resolution = 0.7)
Fang.QC <- RunUMAP(Fang.QC, dims = 1:20)
DimPlot(Fang.QC, reduction = "umap", label = T)
FeaturePlot(Fang.QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
Fang.QC = RenameIdents(Fang.QC, '0' = "Alpha", '1' = "Beta",
                       '2' = "Beta", '3' = "Alpha", '4' = "Beta", 
                       '5' = "Alpha", '6'="Alpha", '7'="CHGA-",
                       '8'= "CHGA-",'9'="Delta", '10'="PPY", '11'="Beta",
                       '12'= "CHGA-",'13' = "Alpha", '14'="Prolif")
Fang.QC$first_identity = Fang.QC@active.ident

## UMAP
DimPlot(Fang.QC, reduction = "umap", label = F, cols = c("#377EB8", "#E41A1C", "grey", "#984EA3", "#1B9E77", "#7570B3"))+NoLegend()

## Subsetting CHGA+ Cells
fang_endocrine = subset(x = Fang.QC, identity = c('Alpha', 'Beta', 'Delta', 'PPY', 'Prolif'))


# Downloading Adult human islet raw data from Xin et al. 2018 (DOI:10.1016/j.cmet.2016.08.018)
# GSE114297 (GSM: GSM3138939, GSM3138940, GSM3138941, GSM3138942, GSM3138943, GSM3138944, GSM3138945, GSM3138946, GSM3138947, GSM3138948, & GSM3138950)

## Download Raw Data, Editing Gene Names, Creating Seurat Object, & Assigning Metadata
setwd("")
Xin1.counts <- Read10X("")
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin1.counts)[colnames(Xin1.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin1.counts)[rownames(Xin1.counts) == "TCEB2"] ="ELOB"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SHFM1"] ="SEM1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "USMG5"] ="ATP5MD"
rownames(Xin1.counts)[rownames(Xin1.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin1.counts)[rownames(Xin1.counts) == "GNB2L1"] ="RACK1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "MRP63"] ="MRPL57"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SF3B14"] ="SF3B6"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin1.counts)[rownames(Xin1.counts) == "CTSL1"] ="CTSL"
rownames(Xin1.counts)[rownames(Xin1.counts) == "TCEB3"] ="ELOA"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C8orf59"] ="RBIS"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C19orf24"] ="FAM174C"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C6orf1"] ="SMIM29"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin1.counts)[rownames(Xin1.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin1.counts)[rownames(Xin1.counts) == "KIAA1715"] ="LNPK"
rownames(Xin1.counts)[rownames(Xin1.counts) == "RTFDC1"] ="RTF2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM134A"] ="RETREG2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "STRA13"] ="CENPX"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C9orf3"] ="AOPEP"
rownames(Xin1.counts)[rownames(Xin1.counts) == "HBT8"] ="PWAR6"
rownames(Xin1.counts)[rownames(Xin1.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin1.counts)[rownames(Xin1.counts) == "WRB"] ="GET1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ASNA1"] ="GET3"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM63B"] ="MINDY2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "WHSC1"] ="NSD2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "APOPT1"] ="COA8"
rownames(Xin1.counts)[rownames(Xin1.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin1.counts)[rownames(Xin1.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "TMEM261"] ="DMAC1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "PPP2R4"] ="PTPA"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C19orf60"] ="REX1BD"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin1.counts)[rownames(Xin1.counts) == "LHFP"] ="LHFPL6"
rownames(Xin1.counts)[rownames(Xin1.counts) == "MESDC2"] ="MESD"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM127B"] ="RTL8A"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C7orf55"] ="FMC1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "LINC00493"] ="SMIM26"
rownames(Xin1.counts)[rownames(Xin1.counts) == "WBSCR22"] ="BUD23"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM105A"] ="OTULINL"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM127A"] ="RTL8C"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C14orf1"] ="ERG28"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "VIMP"] ="SELENOS"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C16orf13"] ="METTL26"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SELM"] ="SELENOM"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "MLLT4"] ="AFDN"
rownames(Xin1.counts)[rownames(Xin1.counts) == "WBP5"] ="TCEAL9"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C19orf70"] ="MICOS13"
rownames(Xin1.counts)[rownames(Xin1.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin1.counts)[rownames(Xin1.counts) == "TCEB1"] ="ELOC"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SEP15"] ="SELENOF"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM96B"] ="CIA02B"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin1.counts)[rownames(Xin1.counts) == "APOA1BP"] ="NAXE"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM46A"] ="TENT5A"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C14orf166"] ="RTRAF"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin1.counts)[rownames(Xin1.counts) == "MINOS1"] ="MICOS10"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C11orf31"] ="SELENOH"
rownames(Xin1.counts)[rownames(Xin1.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin1.counts)[rownames(Xin1.counts) == "LINC01420"] ="NBDY"
rownames(Xin1.counts)[rownames(Xin1.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C7orf73"] ="STMP1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C19orf43"] ="TRIR"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SELT"] ="SELENOT"
rownames(Xin1.counts)[rownames(Xin1.counts) == "AES"] ="TLE5"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SEPW1"] ="SELENOW"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SELK"] ="SELENOK"
rownames(Xin1.counts)[rownames(Xin1.counts) == "HN1"] ="JPT1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin1.counts)[rownames(Xin1.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin1.counts)[rownames(Xin1.counts) == "MYEOV2"] ="COPS9"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C17or76"] ="SNHG29"
rownames(Xin1.counts)[rownames(Xin1.counts) == "MNF1"] ="UQCC2"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C9orf123"] ="DMAC1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C19orf77"] ="SMIM24"
rownames(Xin1.counts)[rownames(Xin1.counts) == "CNIH"] ="CNIH1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin1.counts)[rownames(Xin1.counts) == "NEURL"] ="NEURL1"
rownames(Xin1.counts)[rownames(Xin1.counts) == "LARGE1"] ="LARGE"
rownames(Xin1.counts)[rownames(Xin1.counts) == "C7orf41"] ="MTURN"
Xin1.data <- CreateSeuratObject(counts = Xin1.counts, project = "Adult-X1")
original_dataset = rep("Adult_cadaveric_xin_sample_1", 1744)
age = rep(32,1744)
BMI = rep(24.9,1744)
gender = rep("female", 1744)
Xin1.data$original_dataset = original_dataset
Xin1.data$age = age
Xin1.data$BMI = BMI
Xin1.data$gender = gender

Xin2.counts <- Read10X("")
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin2.counts)[colnames(Xin2.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin2.counts)[rownames(Xin2.counts) == "TCEB2"] ="ELOB"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SHFM1"] ="SEM1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "USMG5"] ="ATP5MD"
rownames(Xin2.counts)[rownames(Xin2.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin2.counts)[rownames(Xin2.counts) == "GNB2L1"] ="RACK1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "MRP63"] ="MRPL57"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SF3B14"] ="SF3B6"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin2.counts)[rownames(Xin2.counts) == "CTSL1"] ="CTSL"
rownames(Xin2.counts)[rownames(Xin2.counts) == "TCEB3"] ="ELOA"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C8orf59"] ="RBIS"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C19orf24"] ="FAM174C"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C6orf1"] ="SMIM29"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin2.counts)[rownames(Xin2.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin2.counts)[rownames(Xin2.counts) == "KIAA1715"] ="LNPK"
rownames(Xin2.counts)[rownames(Xin2.counts) == "RTFDC1"] ="RTF2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM134A"] ="RETREG2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "STRA13"] ="CENPX"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C9orf3"] ="AOPEP"
rownames(Xin2.counts)[rownames(Xin2.counts) == "HBT8"] ="PWAR6"
rownames(Xin2.counts)[rownames(Xin2.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin2.counts)[rownames(Xin2.counts) == "WRB"] ="GET1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ASNA1"] ="GET3"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM63B"] ="MINDY2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "WHSC1"] ="NSD2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "APOPT1"] ="COA8"
rownames(Xin2.counts)[rownames(Xin2.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin2.counts)[rownames(Xin2.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "TMEM261"] ="DMAC1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "PPP2R4"] ="PTPA"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C19orf60"] ="REX1BD"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin2.counts)[rownames(Xin2.counts) == "LHFP"] ="LHFPL6"
rownames(Xin2.counts)[rownames(Xin2.counts) == "MESDC2"] ="MESD"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM127B"] ="RTL8A"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C7orf55"] ="FMC1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "LINC00493"] ="SMIM26"
rownames(Xin2.counts)[rownames(Xin2.counts) == "WBSCR22"] ="BUD23"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM105A"] ="OTULINL"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM127A"] ="RTL8C"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C14orf1"] ="ERG28"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "VIMP"] ="SELENOS"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C16orf13"] ="METTL26"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SELM"] ="SELENOM"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "MLLT4"] ="AFDN"
rownames(Xin2.counts)[rownames(Xin2.counts) == "WBP5"] ="TCEAL9"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C19orf70"] ="MICOS13"
rownames(Xin2.counts)[rownames(Xin2.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin2.counts)[rownames(Xin2.counts) == "TCEB1"] ="ELOC"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SEP15"] ="SELENOF"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM96B"] ="CIA02B"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin2.counts)[rownames(Xin2.counts) == "APOA1BP"] ="NAXE"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM46A"] ="TENT5A"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C14orf166"] ="RTRAF"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin2.counts)[rownames(Xin2.counts) == "MINOS1"] ="MICOS10"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C11orf31"] ="SELENOH"
rownames(Xin2.counts)[rownames(Xin2.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin2.counts)[rownames(Xin2.counts) == "LINC01420"] ="NBDY"
rownames(Xin2.counts)[rownames(Xin2.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C7orf73"] ="STMP1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C19orf43"] ="TRIR"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SELT"] ="SELENOT"
rownames(Xin2.counts)[rownames(Xin2.counts) == "AES"] ="TLE5"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SEPW1"] ="SELENOW"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SELK"] ="SELENOK"
rownames(Xin2.counts)[rownames(Xin2.counts) == "HN1"] ="JPT1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin2.counts)[rownames(Xin2.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin2.counts)[rownames(Xin2.counts) == "MYEOV2"] ="COPS9"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C17or76"] ="SNHG29"
rownames(Xin2.counts)[rownames(Xin2.counts) == "MNF1"] ="UQCC2"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C9orf123"] ="DMAC1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C19orf77"] ="SMIM24"
rownames(Xin2.counts)[rownames(Xin2.counts) == "CNIH"] ="CNIH1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin2.counts)[rownames(Xin2.counts) == "NEURL"] ="NEURL1"
rownames(Xin2.counts)[rownames(Xin2.counts) == "LARGE1"] ="LARGE"
rownames(Xin2.counts)[rownames(Xin2.counts) == "C7orf41"] ="MTURN"
Xin2.data <- CreateSeuratObject(counts = Xin2.counts, project = "Adult-X2")
original_dataset = rep("Adult_cadaveric_xin_sample_2", 933)
age = rep(23,933)
BMI = rep(24.8,933)
gender = rep("male", 933)
Xin2.data$original_dataset = original_dataset
Xin2.data$age = age
Xin2.data$BMI = BMI
Xin2.data$gender = gender

Xin3.counts <- Read10X("")
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin3.counts)[colnames(Xin3.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin3.counts)[rownames(Xin3.counts) == "TCEB2"] ="ELOB"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SHFM1"] ="SEM1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "USMG5"] ="ATP5MD"
rownames(Xin3.counts)[rownames(Xin3.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin3.counts)[rownames(Xin3.counts) == "GNB2L1"] ="RACK1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "MRP63"] ="MRPL57"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SF3B14"] ="SF3B6"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin3.counts)[rownames(Xin3.counts) == "CTSL1"] ="CTSL"
rownames(Xin3.counts)[rownames(Xin3.counts) == "TCEB3"] ="ELOA"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C8orf59"] ="RBIS"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C19orf24"] ="FAM174C"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C6orf1"] ="SMIM29"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin3.counts)[rownames(Xin3.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin3.counts)[rownames(Xin3.counts) == "KIAA1715"] ="LNPK"
rownames(Xin3.counts)[rownames(Xin3.counts) == "RTFDC1"] ="RTF2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM134A"] ="RETREG2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "STRA13"] ="CENPX"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C9orf3"] ="AOPEP"
rownames(Xin3.counts)[rownames(Xin3.counts) == "HBT8"] ="PWAR6"
rownames(Xin3.counts)[rownames(Xin3.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin3.counts)[rownames(Xin3.counts) == "WRB"] ="GET1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ASNA1"] ="GET3"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM63B"] ="MINDY2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "WHSC1"] ="NSD2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "APOPT1"] ="COA8"
rownames(Xin3.counts)[rownames(Xin3.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin3.counts)[rownames(Xin3.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "TMEM261"] ="DMAC1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "PPP2R4"] ="PTPA"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C19orf60"] ="REX1BD"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin3.counts)[rownames(Xin3.counts) == "LHFP"] ="LHFPL6"
rownames(Xin3.counts)[rownames(Xin3.counts) == "MESDC2"] ="MESD"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM127B"] ="RTL8A"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C7orf55"] ="FMC1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "LINC00493"] ="SMIM26"
rownames(Xin3.counts)[rownames(Xin3.counts) == "WBSCR22"] ="BUD23"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM105A"] ="OTULINL"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM127A"] ="RTL8C"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C14orf1"] ="ERG28"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "VIMP"] ="SELENOS"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C16orf13"] ="METTL26"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SELM"] ="SELENOM"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "MLLT4"] ="AFDN"
rownames(Xin3.counts)[rownames(Xin3.counts) == "WBP5"] ="TCEAL9"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C19orf70"] ="MICOS13"
rownames(Xin3.counts)[rownames(Xin3.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin3.counts)[rownames(Xin3.counts) == "TCEB1"] ="ELOC"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SEP15"] ="SELENOF"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM96B"] ="CIA02B"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin3.counts)[rownames(Xin3.counts) == "APOA1BP"] ="NAXE"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM46A"] ="TENT5A"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C14orf166"] ="RTRAF"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin3.counts)[rownames(Xin3.counts) == "MINOS1"] ="MICOS10"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C11orf31"] ="SELENOH"
rownames(Xin3.counts)[rownames(Xin3.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin3.counts)[rownames(Xin3.counts) == "LINC01420"] ="NBDY"
rownames(Xin3.counts)[rownames(Xin3.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C7orf73"] ="STMP1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C19orf43"] ="TRIR"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SELT"] ="SELENOT"
rownames(Xin3.counts)[rownames(Xin3.counts) == "AES"] ="TLE5"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SEPW1"] ="SELENOW"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SELK"] ="SELENOK"
rownames(Xin3.counts)[rownames(Xin3.counts) == "HN1"] ="JPT1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin3.counts)[rownames(Xin3.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin3.counts)[rownames(Xin3.counts) == "MYEOV2"] ="COPS9"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C17or76"] ="SNHG29"
rownames(Xin3.counts)[rownames(Xin3.counts) == "MNF1"] ="UQCC2"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C9orf123"] ="DMAC1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C19orf77"] ="SMIM24"
rownames(Xin3.counts)[rownames(Xin3.counts) == "CNIH"] ="CNIH1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin3.counts)[rownames(Xin3.counts) == "NEURL"] ="NEURL1"
rownames(Xin3.counts)[rownames(Xin3.counts) == "LARGE1"] ="LARGE"
rownames(Xin3.counts)[rownames(Xin3.counts) == "C7orf41"] ="MTURN"
Xin3.data <- CreateSeuratObject(counts = Xin3.counts, project = "Adult-X3")
original_dataset = rep("Adult_cadaveric_xin_sample_3", 2009)
age = rep(52,2009)
BMI = rep(30,2009)
gender = rep("male", 2009)
Xin3.data$original_dataset = original_dataset
Xin3.data$age = age
Xin3.data$BMI = BMI
Xin3.data$gender = gender

Xin4.counts <- Read10X("")
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin4.counts)[colnames(Xin4.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin4.counts)[rownames(Xin4.counts) == "TCEB2"] ="ELOB"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SHFM1"] ="SEM1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "USMG5"] ="ATP5MD"
rownames(Xin4.counts)[rownames(Xin4.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin4.counts)[rownames(Xin4.counts) == "GNB2L1"] ="RACK1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "MRP63"] ="MRPL57"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SF3B14"] ="SF3B6"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin4.counts)[rownames(Xin4.counts) == "CTSL1"] ="CTSL"
rownames(Xin4.counts)[rownames(Xin4.counts) == "TCEB3"] ="ELOA"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C8orf59"] ="RBIS"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C19orf24"] ="FAM174C"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C6orf1"] ="SMIM29"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin4.counts)[rownames(Xin4.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin4.counts)[rownames(Xin4.counts) == "KIAA1715"] ="LNPK"
rownames(Xin4.counts)[rownames(Xin4.counts) == "RTFDC1"] ="RTF2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM134A"] ="RETREG2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "STRA13"] ="CENPX"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C9orf3"] ="AOPEP"
rownames(Xin4.counts)[rownames(Xin4.counts) == "HBT8"] ="PWAR6"
rownames(Xin4.counts)[rownames(Xin4.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin4.counts)[rownames(Xin4.counts) == "WRB"] ="GET1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ASNA1"] ="GET3"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM63B"] ="MINDY2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "WHSC1"] ="NSD2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "APOPT1"] ="COA8"
rownames(Xin4.counts)[rownames(Xin4.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin4.counts)[rownames(Xin4.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "TMEM261"] ="DMAC1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "PPP2R4"] ="PTPA"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C19orf60"] ="REX1BD"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin4.counts)[rownames(Xin4.counts) == "LHFP"] ="LHFPL6"
rownames(Xin4.counts)[rownames(Xin4.counts) == "MESDC2"] ="MESD"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM127B"] ="RTL8A"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C7orf55"] ="FMC1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "LINC00493"] ="SMIM26"
rownames(Xin4.counts)[rownames(Xin4.counts) == "WBSCR22"] ="BUD23"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM105A"] ="OTULINL"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM127A"] ="RTL8C"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C14orf1"] ="ERG28"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "VIMP"] ="SELENOS"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C16orf13"] ="METTL26"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SELM"] ="SELENOM"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "MLLT4"] ="AFDN"
rownames(Xin4.counts)[rownames(Xin4.counts) == "WBP5"] ="TCEAL9"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C19orf70"] ="MICOS13"
rownames(Xin4.counts)[rownames(Xin4.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin4.counts)[rownames(Xin4.counts) == "TCEB1"] ="ELOC"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SEP15"] ="SELENOF"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM96B"] ="CIA02B"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin4.counts)[rownames(Xin4.counts) == "APOA1BP"] ="NAXE"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM46A"] ="TENT5A"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C14orf166"] ="RTRAF"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin4.counts)[rownames(Xin4.counts) == "MINOS1"] ="MICOS10"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C11orf31"] ="SELENOH"
rownames(Xin4.counts)[rownames(Xin4.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin4.counts)[rownames(Xin4.counts) == "LINC01420"] ="NBDY"
rownames(Xin4.counts)[rownames(Xin4.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C7orf73"] ="STMP1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C19orf43"] ="TRIR"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SELT"] ="SELENOT"
rownames(Xin4.counts)[rownames(Xin4.counts) == "AES"] ="TLE5"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SEPW1"] ="SELENOW"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SELK"] ="SELENOK"
rownames(Xin4.counts)[rownames(Xin4.counts) == "HN1"] ="JPT1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin4.counts)[rownames(Xin4.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin4.counts)[rownames(Xin4.counts) == "MYEOV2"] ="COPS9"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C17or76"] ="SNHG29"
rownames(Xin4.counts)[rownames(Xin4.counts) == "MNF1"] ="UQCC2"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C9orf123"] ="DMAC1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C19orf77"] ="SMIM24"
rownames(Xin4.counts)[rownames(Xin4.counts) == "CNIH"] ="CNIH1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin4.counts)[rownames(Xin4.counts) == "NEURL"] ="NEURL1"
rownames(Xin4.counts)[rownames(Xin4.counts) == "LARGE1"] ="LARGE"
rownames(Xin4.counts)[rownames(Xin4.counts) == "C7orf41"] ="MTURN"
Xin4.data <- CreateSeuratObject(counts = Xin4.counts, project = "Adult-X4")
original_dataset = rep("Adult_cadaveric_xin_sample_4", 1534)
age = rep(52,1534)
BMI = rep(22,1534)
gender = rep("male", 1534)
Xin4.data$original_dataset = original_dataset
Xin4.data$age = age
Xin4.data$BMI = BMI
Xin4.data$gender = gender

Xin5.counts <- Read10X("")
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin5.counts)[colnames(Xin5.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin5.counts)[rownames(Xin5.counts) == "TCEB2"] ="ELOB"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SHFM1"] ="SEM1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "USMG5"] ="ATP5MD"
rownames(Xin5.counts)[rownames(Xin5.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin5.counts)[rownames(Xin5.counts) == "GNB2L1"] ="RACK1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "MRP63"] ="MRPL57"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SF3B14"] ="SF3B6"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin5.counts)[rownames(Xin5.counts) == "CTSL1"] ="CTSL"
rownames(Xin5.counts)[rownames(Xin5.counts) == "TCEB3"] ="ELOA"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C8orf59"] ="RBIS"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C19orf24"] ="FAM174C"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C6orf1"] ="SMIM29"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin5.counts)[rownames(Xin5.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin5.counts)[rownames(Xin5.counts) == "KIAA1715"] ="LNPK"
rownames(Xin5.counts)[rownames(Xin5.counts) == "RTFDC1"] ="RTF2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM134A"] ="RETREG2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "STRA13"] ="CENPX"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C9orf3"] ="AOPEP"
rownames(Xin5.counts)[rownames(Xin5.counts) == "HBT8"] ="PWAR6"
rownames(Xin5.counts)[rownames(Xin5.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin5.counts)[rownames(Xin5.counts) == "WRB"] ="GET1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ASNA1"] ="GET3"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM63B"] ="MINDY2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "WHSC1"] ="NSD2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "APOPT1"] ="COA8"
rownames(Xin5.counts)[rownames(Xin5.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin5.counts)[rownames(Xin5.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "TMEM261"] ="DMAC1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "PPP2R4"] ="PTPA"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C19orf60"] ="REX1BD"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin5.counts)[rownames(Xin5.counts) == "LHFP"] ="LHFPL6"
rownames(Xin5.counts)[rownames(Xin5.counts) == "MESDC2"] ="MESD"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM127B"] ="RTL8A"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C7orf55"] ="FMC1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "LINC00493"] ="SMIM26"
rownames(Xin5.counts)[rownames(Xin5.counts) == "WBSCR22"] ="BUD23"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM105A"] ="OTULINL"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM127A"] ="RTL8C"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C14orf1"] ="ERG28"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "VIMP"] ="SELENOS"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C16orf13"] ="METTL26"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SELM"] ="SELENOM"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "MLLT4"] ="AFDN"
rownames(Xin5.counts)[rownames(Xin5.counts) == "WBP5"] ="TCEAL9"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C19orf70"] ="MICOS13"
rownames(Xin5.counts)[rownames(Xin5.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin5.counts)[rownames(Xin5.counts) == "TCEB1"] ="ELOC"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SEP15"] ="SELENOF"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM96B"] ="CIA02B"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin5.counts)[rownames(Xin5.counts) == "APOA1BP"] ="NAXE"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM46A"] ="TENT5A"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C14orf166"] ="RTRAF"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin5.counts)[rownames(Xin5.counts) == "MINOS1"] ="MICOS10"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C11orf31"] ="SELENOH"
rownames(Xin5.counts)[rownames(Xin5.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin5.counts)[rownames(Xin5.counts) == "LINC01420"] ="NBDY"
rownames(Xin5.counts)[rownames(Xin5.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C7orf73"] ="STMP1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C19orf43"] ="TRIR"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SELT"] ="SELENOT"
rownames(Xin5.counts)[rownames(Xin5.counts) == "AES"] ="TLE5"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SEPW1"] ="SELENOW"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SELK"] ="SELENOK"
rownames(Xin5.counts)[rownames(Xin5.counts) == "HN1"] ="JPT1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin5.counts)[rownames(Xin5.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin5.counts)[rownames(Xin5.counts) == "MYEOV2"] ="COPS9"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C17or76"] ="SNHG29"
rownames(Xin5.counts)[rownames(Xin5.counts) == "MNF1"] ="UQCC2"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C9orf123"] ="DMAC1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C19orf77"] ="SMIM24"
rownames(Xin5.counts)[rownames(Xin5.counts) == "CNIH"] ="CNIH1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin5.counts)[rownames(Xin5.counts) == "NEURL"] ="NEURL1"
rownames(Xin5.counts)[rownames(Xin5.counts) == "LARGE1"] ="LARGE"
rownames(Xin5.counts)[rownames(Xin5.counts) == "C7orf41"] ="MTURN"
Xin5.data <- CreateSeuratObject(counts = Xin5.counts, project = "Adult-X5")
original_dataset = rep("Adult_cadaveric_xin_sample_5", 1068)
age = rep(49,1068)
BMI = rep(26.5,1068)
gender = rep("male", 1068)
Xin5.data$original_dataset = original_dataset
Xin5.data$age = age
Xin5.data$BMI = BMI
Xin5.data$gender = gender

Xin6.counts <- Read10X("")
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin6.counts)[colnames(Xin6.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin6.counts)[rownames(Xin6.counts) == "TCEB2"] ="ELOB"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SHFM1"] ="SEM1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "USMG5"] ="ATP5MD"
rownames(Xin6.counts)[rownames(Xin6.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin6.counts)[rownames(Xin6.counts) == "GNB2L1"] ="RACK1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "MRP63"] ="MRPL57"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SF3B14"] ="SF3B6"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin6.counts)[rownames(Xin6.counts) == "CTSL1"] ="CTSL"
rownames(Xin6.counts)[rownames(Xin6.counts) == "TCEB3"] ="ELOA"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C8orf59"] ="RBIS"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C19orf24"] ="FAM174C"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C6orf1"] ="SMIM29"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin6.counts)[rownames(Xin6.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin6.counts)[rownames(Xin6.counts) == "KIAA1715"] ="LNPK"
rownames(Xin6.counts)[rownames(Xin6.counts) == "RTFDC1"] ="RTF2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM134A"] ="RETREG2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "STRA13"] ="CENPX"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C9orf3"] ="AOPEP"
rownames(Xin6.counts)[rownames(Xin6.counts) == "HBT8"] ="PWAR6"
rownames(Xin6.counts)[rownames(Xin6.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin6.counts)[rownames(Xin6.counts) == "WRB"] ="GET1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ASNA1"] ="GET3"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM63B"] ="MINDY2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "WHSC1"] ="NSD2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "APOPT1"] ="COA8"
rownames(Xin6.counts)[rownames(Xin6.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin6.counts)[rownames(Xin6.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "TMEM261"] ="DMAC1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "PPP2R4"] ="PTPA"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C19orf60"] ="REX1BD"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin6.counts)[rownames(Xin6.counts) == "LHFP"] ="LHFPL6"
rownames(Xin6.counts)[rownames(Xin6.counts) == "MESDC2"] ="MESD"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM127B"] ="RTL8A"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C7orf55"] ="FMC1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "LINC00493"] ="SMIM26"
rownames(Xin6.counts)[rownames(Xin6.counts) == "WBSCR22"] ="BUD23"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM105A"] ="OTULINL"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM127A"] ="RTL8C"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C14orf1"] ="ERG28"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "VIMP"] ="SELENOS"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C16orf13"] ="METTL26"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SELM"] ="SELENOM"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "MLLT4"] ="AFDN"
rownames(Xin6.counts)[rownames(Xin6.counts) == "WBP5"] ="TCEAL9"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C19orf70"] ="MICOS13"
rownames(Xin6.counts)[rownames(Xin6.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin6.counts)[rownames(Xin6.counts) == "TCEB1"] ="ELOC"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SEP15"] ="SELENOF"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM96B"] ="CIA02B"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin6.counts)[rownames(Xin6.counts) == "APOA1BP"] ="NAXE"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM46A"] ="TENT5A"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C14orf166"] ="RTRAF"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin6.counts)[rownames(Xin6.counts) == "MINOS1"] ="MICOS10"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C11orf31"] ="SELENOH"
rownames(Xin6.counts)[rownames(Xin6.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin6.counts)[rownames(Xin6.counts) == "LINC01420"] ="NBDY"
rownames(Xin6.counts)[rownames(Xin6.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C7orf73"] ="STMP1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C19orf43"] ="TRIR"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SELT"] ="SELENOT"
rownames(Xin6.counts)[rownames(Xin6.counts) == "AES"] ="TLE5"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SEPW1"] ="SELENOW"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SELK"] ="SELENOK"
rownames(Xin6.counts)[rownames(Xin6.counts) == "HN1"] ="JPT1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin6.counts)[rownames(Xin6.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin6.counts)[rownames(Xin6.counts) == "MYEOV2"] ="COPS9"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C17or76"] ="SNHG29"
rownames(Xin6.counts)[rownames(Xin6.counts) == "MNF1"] ="UQCC2"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C9orf123"] ="DMAC1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C19orf77"] ="SMIM24"
rownames(Xin6.counts)[rownames(Xin6.counts) == "CNIH"] ="CNIH1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin6.counts)[rownames(Xin6.counts) == "NEURL"] ="NEURL1"
rownames(Xin6.counts)[rownames(Xin6.counts) == "LARGE1"] ="LARGE"
rownames(Xin6.counts)[rownames(Xin6.counts) == "C7orf41"] ="MTURN"
Xin6.data <- CreateSeuratObject(counts = Xin6.counts, project = "Adult-X6")
original_dataset = rep("Adult_cadaveric_xin_sample_6", 410)
age = rep(37,410)
BMI = rep(22.9,410)
gender = rep("female", 410)
Xin6.data$original_dataset = original_dataset
Xin6.data$age = age
Xin6.data$BMI = BMI
Xin6.data$gender = gender

Xin7.counts <- Read10X("")
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin7.counts)[colnames(Xin7.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin7.counts)[rownames(Xin7.counts) == "TCEB2"] ="ELOB"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SHFM1"] ="SEM1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "USMG5"] ="ATP5MD"
rownames(Xin7.counts)[rownames(Xin7.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin7.counts)[rownames(Xin7.counts) == "GNB2L1"] ="RACK1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "MRP63"] ="MRPL57"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SF3B14"] ="SF3B6"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin7.counts)[rownames(Xin7.counts) == "CTSL1"] ="CTSL"
rownames(Xin7.counts)[rownames(Xin7.counts) == "TCEB3"] ="ELOA"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C8orf59"] ="RBIS"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C19orf24"] ="FAM174C"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C6orf1"] ="SMIM29"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin7.counts)[rownames(Xin7.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin7.counts)[rownames(Xin7.counts) == "KIAA1715"] ="LNPK"
rownames(Xin7.counts)[rownames(Xin7.counts) == "RTFDC1"] ="RTF2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM134A"] ="RETREG2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "STRA13"] ="CENPX"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C9orf3"] ="AOPEP"
rownames(Xin7.counts)[rownames(Xin7.counts) == "HBT8"] ="PWAR6"
rownames(Xin7.counts)[rownames(Xin7.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin7.counts)[rownames(Xin7.counts) == "WRB"] ="GET1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ASNA1"] ="GET3"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM63B"] ="MINDY2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "WHSC1"] ="NSD2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "APOPT1"] ="COA8"
rownames(Xin7.counts)[rownames(Xin7.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin7.counts)[rownames(Xin7.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "TMEM261"] ="DMAC1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "PPP2R4"] ="PTPA"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C19orf60"] ="REX1BD"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin7.counts)[rownames(Xin7.counts) == "LHFP"] ="LHFPL6"
rownames(Xin7.counts)[rownames(Xin7.counts) == "MESDC2"] ="MESD"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM127B"] ="RTL8A"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C7orf55"] ="FMC1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "LINC00493"] ="SMIM26"
rownames(Xin7.counts)[rownames(Xin7.counts) == "WBSCR22"] ="BUD23"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM105A"] ="OTULINL"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM127A"] ="RTL8C"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C14orf1"] ="ERG28"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "VIMP"] ="SELENOS"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C16orf13"] ="METTL26"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SELM"] ="SELENOM"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "MLLT4"] ="AFDN"
rownames(Xin7.counts)[rownames(Xin7.counts) == "WBP5"] ="TCEAL9"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C19orf70"] ="MICOS13"
rownames(Xin7.counts)[rownames(Xin7.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin7.counts)[rownames(Xin7.counts) == "TCEB1"] ="ELOC"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SEP15"] ="SELENOF"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM96B"] ="CIA02B"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin7.counts)[rownames(Xin7.counts) == "APOA1BP"] ="NAXE"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM46A"] ="TENT5A"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C14orf166"] ="RTRAF"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin7.counts)[rownames(Xin7.counts) == "MINOS1"] ="MICOS10"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C11orf31"] ="SELENOH"
rownames(Xin7.counts)[rownames(Xin7.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin7.counts)[rownames(Xin7.counts) == "LINC01420"] ="NBDY"
rownames(Xin7.counts)[rownames(Xin7.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C7orf73"] ="STMP1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C19orf43"] ="TRIR"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SELT"] ="SELENOT"
rownames(Xin7.counts)[rownames(Xin7.counts) == "AES"] ="TLE5"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SEPW1"] ="SELENOW"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SELK"] ="SELENOK"
rownames(Xin7.counts)[rownames(Xin7.counts) == "HN1"] ="JPT1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin7.counts)[rownames(Xin7.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin7.counts)[rownames(Xin7.counts) == "MYEOV2"] ="COPS9"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C17or76"] ="SNHG29"
rownames(Xin7.counts)[rownames(Xin7.counts) == "MNF1"] ="UQCC2"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C9orf123"] ="DMAC1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C19orf77"] ="SMIM24"
rownames(Xin7.counts)[rownames(Xin7.counts) == "CNIH"] ="CNIH1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin7.counts)[rownames(Xin7.counts) == "NEURL"] ="NEURL1"
rownames(Xin7.counts)[rownames(Xin7.counts) == "LARGE1"] ="LARGE"
rownames(Xin7.counts)[rownames(Xin7.counts) == "C7orf41"] ="MTURN"
Xin7.data <- CreateSeuratObject(counts = Xin7.counts, project = "Adult-X7")
original_dataset = rep("Adult_cadaveric_xin_sample_7", 1883)
age = rep(39,1883)
BMI = rep(23.6,1883)
gender = rep("male", 1883)
Xin7.data$original_dataset = original_dataset
Xin7.data$age = age
Xin7.data$BMI = BMI
Xin7.data$gender = gender

Xin8.counts <- Read10X("")
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin8.counts)[colnames(Xin8.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin8.counts)[rownames(Xin8.counts) == "TCEB2"] ="ELOB"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SHFM1"] ="SEM1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "USMG5"] ="ATP5MD"
rownames(Xin8.counts)[rownames(Xin8.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin8.counts)[rownames(Xin8.counts) == "GNB2L1"] ="RACK1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "MRP63"] ="MRPL57"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SF3B14"] ="SF3B6"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin8.counts)[rownames(Xin8.counts) == "CTSL1"] ="CTSL"
rownames(Xin8.counts)[rownames(Xin8.counts) == "TCEB3"] ="ELOA"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C8orf59"] ="RBIS"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C19orf24"] ="FAM174C"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C6orf1"] ="SMIM29"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin8.counts)[rownames(Xin8.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin8.counts)[rownames(Xin8.counts) == "KIAA1715"] ="LNPK"
rownames(Xin8.counts)[rownames(Xin8.counts) == "RTFDC1"] ="RTF2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM134A"] ="RETREG2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "STRA13"] ="CENPX"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C9orf3"] ="AOPEP"
rownames(Xin8.counts)[rownames(Xin8.counts) == "HBT8"] ="PWAR6"
rownames(Xin8.counts)[rownames(Xin8.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin8.counts)[rownames(Xin8.counts) == "WRB"] ="GET1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ASNA1"] ="GET3"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM63B"] ="MINDY2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "WHSC1"] ="NSD2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "APOPT1"] ="COA8"
rownames(Xin8.counts)[rownames(Xin8.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin8.counts)[rownames(Xin8.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "TMEM261"] ="DMAC1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "PPP2R4"] ="PTPA"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C19orf60"] ="REX1BD"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin8.counts)[rownames(Xin8.counts) == "LHFP"] ="LHFPL6"
rownames(Xin8.counts)[rownames(Xin8.counts) == "MESDC2"] ="MESD"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM127B"] ="RTL8A"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C7orf55"] ="FMC1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "LINC00493"] ="SMIM26"
rownames(Xin8.counts)[rownames(Xin8.counts) == "WBSCR22"] ="BUD23"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM105A"] ="OTULINL"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM127A"] ="RTL8C"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C14orf1"] ="ERG28"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "VIMP"] ="SELENOS"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C16orf13"] ="METTL26"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SELM"] ="SELENOM"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "MLLT4"] ="AFDN"
rownames(Xin8.counts)[rownames(Xin8.counts) == "WBP5"] ="TCEAL9"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C19orf70"] ="MICOS13"
rownames(Xin8.counts)[rownames(Xin8.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin8.counts)[rownames(Xin8.counts) == "TCEB1"] ="ELOC"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SEP15"] ="SELENOF"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM96B"] ="CIA02B"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin8.counts)[rownames(Xin8.counts) == "APOA1BP"] ="NAXE"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM46A"] ="TENT5A"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C14orf166"] ="RTRAF"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin8.counts)[rownames(Xin8.counts) == "MINOS1"] ="MICOS10"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C11orf31"] ="SELENOH"
rownames(Xin8.counts)[rownames(Xin8.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin8.counts)[rownames(Xin8.counts) == "LINC01420"] ="NBDY"
rownames(Xin8.counts)[rownames(Xin8.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C7orf73"] ="STMP1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C19orf43"] ="TRIR"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SELT"] ="SELENOT"
rownames(Xin8.counts)[rownames(Xin8.counts) == "AES"] ="TLE5"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SEPW1"] ="SELENOW"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SELK"] ="SELENOK"
rownames(Xin8.counts)[rownames(Xin8.counts) == "HN1"] ="JPT1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin8.counts)[rownames(Xin8.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin8.counts)[rownames(Xin8.counts) == "MYEOV2"] ="COPS9"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C17or76"] ="SNHG29"
rownames(Xin8.counts)[rownames(Xin8.counts) == "MNF1"] ="UQCC2"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C9orf123"] ="DMAC1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C19orf77"] ="SMIM24"
rownames(Xin8.counts)[rownames(Xin8.counts) == "CNIH"] ="CNIH1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin8.counts)[rownames(Xin8.counts) == "NEURL"] ="NEURL1"
rownames(Xin8.counts)[rownames(Xin8.counts) == "LARGE1"] ="LARGE"
rownames(Xin8.counts)[rownames(Xin8.counts) == "C7orf41"] ="MTURN"
Xin8.data <- CreateSeuratObject(counts = Xin8.counts, project = "Adult-X8")
original_dataset = rep("Adult_cadaveric_xin_sample_8", 1571)
age = rep(28,1571)
BMI = rep(30.8,1571)
gender = rep("male", 1571)
Xin8.data$original_dataset = original_dataset
Xin8.data$age = age
Xin8.data$BMI = BMI
Xin8.data$gender = gender

Xin9.counts <- Read10X("")
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin9.counts)[colnames(Xin9.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin9.counts)[rownames(Xin9.counts) == "TCEB2"] ="ELOB"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SHFM1"] ="SEM1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "USMG5"] ="ATP5MD"
rownames(Xin9.counts)[rownames(Xin9.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin9.counts)[rownames(Xin9.counts) == "GNB2L1"] ="RACK1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "MRP63"] ="MRPL57"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SF3B14"] ="SF3B6"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin9.counts)[rownames(Xin9.counts) == "CTSL1"] ="CTSL"
rownames(Xin9.counts)[rownames(Xin9.counts) == "TCEB3"] ="ELOA"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C8orf59"] ="RBIS"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C19orf24"] ="FAM174C"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C6orf1"] ="SMIM29"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin9.counts)[rownames(Xin9.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin9.counts)[rownames(Xin9.counts) == "KIAA1715"] ="LNPK"
rownames(Xin9.counts)[rownames(Xin9.counts) == "RTFDC1"] ="RTF2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM134A"] ="RETREG2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "STRA13"] ="CENPX"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C9orf3"] ="AOPEP"
rownames(Xin9.counts)[rownames(Xin9.counts) == "HBT8"] ="PWAR6"
rownames(Xin9.counts)[rownames(Xin9.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin9.counts)[rownames(Xin9.counts) == "WRB"] ="GET1"
rownames(Xin9.counts)[rownames(Xin8.counts) == "ASNA1"] ="GET3"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM63B"] ="MINDY2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "WHSC1"] ="NSD2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "APOPT1"] ="COA8"
rownames(Xin9.counts)[rownames(Xin9.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin9.counts)[rownames(Xin9.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "TMEM261"] ="DMAC1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "PPP2R4"] ="PTPA"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C19orf60"] ="REX1BD"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin9.counts)[rownames(Xin9.counts) == "LHFP"] ="LHFPL6"
rownames(Xin9.counts)[rownames(Xin9.counts) == "MESDC2"] ="MESD"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM127B"] ="RTL8A"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C7orf55"] ="FMC1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "LINC00493"] ="SMIM26"
rownames(Xin9.counts)[rownames(Xin9.counts) == "WBSCR22"] ="BUD23"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM105A"] ="OTULINL"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM127A"] ="RTL8C"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C14orf1"] ="ERG28"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "VIMP"] ="SELENOS"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C16orf13"] ="METTL26"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SELM"] ="SELENOM"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "MLLT4"] ="AFDN"
rownames(Xin9.counts)[rownames(Xin9.counts) == "WBP5"] ="TCEAL9"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C19orf70"] ="MICOS13"
rownames(Xin9.counts)[rownames(Xin9.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin9.counts)[rownames(Xin9.counts) == "TCEB1"] ="ELOC"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SEP15"] ="SELENOF"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM96B"] ="CIA02B"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin9.counts)[rownames(Xin9.counts) == "APOA1BP"] ="NAXE"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM46A"] ="TENT5A"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C14orf166"] ="RTRAF"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin9.counts)[rownames(Xin9.counts) == "MINOS1"] ="MICOS10"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C11orf31"] ="SELENOH"
rownames(Xin9.counts)[rownames(Xin9.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin9.counts)[rownames(Xin9.counts) == "LINC01420"] ="NBDY"
rownames(Xin9.counts)[rownames(Xin9.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C7orf73"] ="STMP1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C19orf43"] ="TRIR"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SELT"] ="SELENOT"
rownames(Xin9.counts)[rownames(Xin9.counts) == "AES"] ="TLE5"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SEPW1"] ="SELENOW"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SELK"] ="SELENOK"
rownames(Xin9.counts)[rownames(Xin9.counts) == "HN1"] ="JPT1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin9.counts)[rownames(Xin9.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin9.counts)[rownames(Xin9.counts) == "MYEOV2"] ="COPS9"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C17or76"] ="SNHG29"
rownames(Xin9.counts)[rownames(Xin9.counts) == "MNF1"] ="UQCC2"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C9orf123"] ="DMAC1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C19orf77"] ="SMIM24"
rownames(Xin9.counts)[rownames(Xin9.counts) == "CNIH"] ="CNIH1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin9.counts)[rownames(Xin9.counts) == "NEURL"] ="NEURL1"
rownames(Xin9.counts)[rownames(Xin9.counts) == "LARGE1"] ="LARGE"
rownames(Xin9.counts)[rownames(Xin9.counts) == "C7orf41"] ="MTURN"
Xin9.data <- CreateSeuratObject(counts = Xin9.counts, project = "Adult-X9")
original_dataset = rep("Adult_cadaveric_xin_sample_9", 2284)
age = rep(49,2284)
BMI = rep(26.2,2284)
gender = rep("male", 2284)
Xin9.data$original_dataset = original_dataset
Xin9.data$age = age
Xin9.data$BMI = BMI
Xin9.data$gender = gender

Xin10.counts <- Read10X("")
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin10.counts)[colnames(Xin10.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin10.counts)[rownames(Xin10.counts) == "TCEB2"] ="ELOB"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SHFM1"] ="SEM1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "USMG5"] ="ATP5MD"
rownames(Xin10.counts)[rownames(Xin10.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin10.counts)[rownames(Xin10.counts) == "GNB2L1"] ="RACK1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "MRP63"] ="MRPL57"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SF3B14"] ="SF3B6"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin10.counts)[rownames(Xin10.counts) == "CTSL1"] ="CTSL"
rownames(Xin10.counts)[rownames(Xin10.counts) == "TCEB3"] ="ELOA"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C8orf59"] ="RBIS"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C19orf24"] ="FAM174C"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C6orf1"] ="SMIM29"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin10.counts)[rownames(Xin10.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin10.counts)[rownames(Xin10.counts) == "KIAA1715"] ="LNPK"
rownames(Xin10.counts)[rownames(Xin10.counts) == "RTFDC1"] ="RTF2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM134A"] ="RETREG2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "STRA13"] ="CENPX"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C9orf3"] ="AOPEP"
rownames(Xin10.counts)[rownames(Xin10.counts) == "HBT8"] ="PWAR6"
rownames(Xin10.counts)[rownames(Xin10.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin10.counts)[rownames(Xin10.counts) == "WRB"] ="GET1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ASNA1"] ="GET3"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM63B"] ="MINDY2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "WHSC1"] ="NSD2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "APOPT1"] ="COA8"
rownames(Xin10.counts)[rownames(Xin10.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin10.counts)[rownames(Xin10.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "TMEM261"] ="DMAC1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "PPP2R4"] ="PTPA"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C19orf60"] ="REX1BD"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin10.counts)[rownames(Xin10.counts) == "LHFP"] ="LHFPL6"
rownames(Xin10.counts)[rownames(Xin10.counts) == "MESDC2"] ="MESD"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM127B"] ="RTL8A"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C7orf55"] ="FMC1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "LINC00493"] ="SMIM26"
rownames(Xin10.counts)[rownames(Xin10.counts) == "WBSCR22"] ="BUD23"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM105A"] ="OTULINL"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM127A"] ="RTL8C"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C14orf1"] ="ERG28"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "VIMP"] ="SELENOS"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C16orf13"] ="METTL26"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SELM"] ="SELENOM"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "MLLT4"] ="AFDN"
rownames(Xin10.counts)[rownames(Xin10.counts) == "WBP5"] ="TCEAL9"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C19orf70"] ="MICOS13"
rownames(Xin10.counts)[rownames(Xin10.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin10.counts)[rownames(Xin10.counts) == "TCEB1"] ="ELOC"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SEP15"] ="SELENOF"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM96B"] ="CIA02B"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin10.counts)[rownames(Xin10.counts) == "APOA1BP"] ="NAXE"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM46A"] ="TENT5A"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C14orf166"] ="RTRAF"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin10.counts)[rownames(Xin10.counts) == "MINOS1"] ="MICOS10"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C11orf31"] ="SELENOH"
rownames(Xin10.counts)[rownames(Xin10.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin10.counts)[rownames(Xin10.counts) == "LINC01420"] ="NBDY"
rownames(Xin10.counts)[rownames(Xin10.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C7orf73"] ="STMP1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C19orf43"] ="TRIR"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SELT"] ="SELENOT"
rownames(Xin10.counts)[rownames(Xin10.counts) == "AES"] ="TLE5"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SEPW1"] ="SELENOW"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SELK"] ="SELENOK"
rownames(Xin10.counts)[rownames(Xin10.counts) == "HN1"] ="JPT1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin10.counts)[rownames(Xin10.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin10.counts)[rownames(Xin10.counts) == "MYEOV2"] ="COPS9"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C17or76"] ="SNHG29"
rownames(Xin10.counts)[rownames(Xin10.counts) == "MNF1"] ="UQCC2"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C9orf123"] ="DMAC1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C19orf77"] ="SMIM24"
rownames(Xin10.counts)[rownames(Xin10.counts) == "CNIH"] ="CNIH1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin10.counts)[rownames(Xin10.counts) == "NEURL"] ="NEURL1"
rownames(Xin10.counts)[rownames(Xin10.counts) == "LARGE1"] ="LARGE"
rownames(Xin10.counts)[rownames(Xin10.counts) == "C7orf41"] ="MTURN"
Xin10.data <- CreateSeuratObject(counts = Xin10.counts, project = "Adult-X10")
original_dataset = rep("Adult_cadaveric_xin_sample_10", 2351)
age = rep(55,2351)
BMI = rep(23.5,2351)
gender = rep("female", 2351)
Xin10.data$original_dataset = original_dataset
Xin10.data$age = age
Xin10.data$BMI = BMI
Xin10.data$gender = gender

Xin12.counts <- Read10X("")
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5I"] ="ATP5ME"
rownames(Xin12.counts)[colnames(Xin12.counts) == "ATP5G2"] ="ATP5MC2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5E"] ="ATP5F1E"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5L"] ="ATP5MG"
rownames(Xin12.counts)[rownames(Xin12.counts) == "TCEB2"] ="ELOB"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SHFM1"] ="SEM1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "USMG5"] ="ATP5MD"
rownames(Xin12.counts)[rownames(Xin12.counts) == "NGFRAP1"] ="BEX3"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM213A"] ="PRXL2A"
rownames(Xin12.counts)[rownames(Xin12.counts) == "GNB2L1"] ="RACK1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "MRP63"] ="MRPL57"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SF3B14"] ="SF3B6"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ZMYM6NB"] ="TMEM35B"
rownames(Xin12.counts)[rownames(Xin12.counts) == "CTSL1"] ="CTSL"
rownames(Xin12.counts)[rownames(Xin12.counts) == "TCEB3"] ="ELOA"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C8orf59"] ="RBIS"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C19orf24"] ="FAM174C"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C20orf24"] ="RAB5IF"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM173A"] ="ANTKMT"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C6orf1"] ="SMIM29"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM57B"] ="TLCD3B"
rownames(Xin12.counts)[rownames(Xin12.counts) == "PVRL2"] ="NECTIN2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM96A"] ="CIAO2A"
rownames(Xin12.counts)[rownames(Xin12.counts) == "KIAA1715"] ="LNPK"
rownames(Xin12.counts)[rownames(Xin12.counts) == "RTFDC1"] ="RTF2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM134A"] ="RETREG2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "STRA13"] ="CENPX"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C9orf3"] ="AOPEP"
rownames(Xin12.counts)[rownames(Xin12.counts) == "HBT8"] ="PWAR6"
rownames(Xin12.counts)[rownames(Xin12.counts) == "PVRL3"] ="NECTIN3"
rownames(Xin12.counts)[rownames(Xin12.counts) == "WRB"] ="GET1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ASNA1"] ="GET3"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM60A"] ="SINHCAF"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM63B"] ="MINDY2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "COL4A3BP"] ="CERT1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "WHSC1"] ="NSD2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "APOPT1"] ="COA8"
rownames(Xin12.counts)[rownames(Xin12.counts) == "TOMM70A"] ="TOMM70"
rownames(Xin12.counts)[rownames(Xin12.counts) == "PVRL1"] ="NECTIN1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "TMEM261"] ="DMAC1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM195A"] ="MCRIP2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "PPP2R4"] ="PTPA"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C11orf73"] ="HIKESHI"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C19orf60"] ="REX1BD"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ZCCHC11"] ="TUT4"
rownames(Xin12.counts)[rownames(Xin12.counts) == "LHFP"] ="LHFPL6"
rownames(Xin12.counts)[rownames(Xin12.counts) == "MESDC2"] ="MESD"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM127B"] ="RTL8A"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C7orf55"] ="FMC1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "LINC00493"] ="SMIM26"
rownames(Xin12.counts)[rownames(Xin12.counts) == "WBSCR22"] ="BUD23"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM105A"] ="OTULINL"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM127A"] ="RTL8C"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C14orf1"] ="ERG28"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SEPT2"] ="SEPTIN2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "VIMP"] ="SELENOS"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C16orf13"] ="METTL26"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SELM"] ="SELENOM"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM195B"] ="MCRIP1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "MLLT4"] ="AFDN"
rownames(Xin12.counts)[rownames(Xin12.counts) == "WBP5"] ="TCEAL9"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SEPT11"] ="SEPTIN11"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C19orf70"] ="MICOS13"
rownames(Xin12.counts)[rownames(Xin12.counts) == "WHSC1L1"] ="NSD3"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5H"] ="ATP5PD"
rownames(Xin12.counts)[rownames(Xin12.counts) == "TCEB1"] ="ELOC"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5C1"] ="ATP5F1C"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SEP15"] ="SELENOF"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM96B"] ="CIA02B"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C17orf89"] ="NDUFAF8"
rownames(Xin12.counts)[rownames(Xin12.counts) == "APOA1BP"] ="NAXE"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5F1"] ="ATP5PB"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM46A"] ="TENT5A"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C14orf166"] ="RTRAF"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SEPT7"] ="SEPTIN7"
rownames(Xin12.counts)[rownames(Xin12.counts) == "MINOS1"] ="MICOS10"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C11orf31"] ="SELENOH"
rownames(Xin12.counts)[rownames(Xin12.counts) == "GLTSCR2"] ="NOP53"
rownames(Xin12.counts)[rownames(Xin12.counts) == "LINC01420"] ="NBDY"
rownames(Xin12.counts)[rownames(Xin12.counts) == "FAM159B"] ="SHISAL2B"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C7orf73"] ="STMP1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C19orf43"] ="TRIR"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SELT"] ="SELENOT"
rownames(Xin12.counts)[rownames(Xin12.counts) == "AES"] ="TLE5"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5G1"] ="ATP5MC1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5B"] ="ATP5F1B"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SEPW1"] ="SELENOW"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5D"] ="ATP5F1D"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5O"] ="ATP5PO"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SELK"] ="SELENOK"
rownames(Xin12.counts)[rownames(Xin12.counts) == "HN1"] ="JPT1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATPIF1"] ="ATP5IF1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5J"] ="ATP5PF"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5A1"] ="ATP5F1A"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5G3"] ="ATP5MC3"
rownames(Xin12.counts)[rownames(Xin12.counts) == "ATP5J2"] ="ATP5MF"
rownames(Xin12.counts)[rownames(Xin12.counts) == "MYEOV2"] ="COPS9"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C14orf2"] ="ATP5MPL"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C17or76"] ="SNHG29"
rownames(Xin12.counts)[rownames(Xin12.counts) == "MNF1"] ="UQCC2"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C9orf123"] ="DMAC1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C19orf77"] ="SMIM24"
rownames(Xin12.counts)[rownames(Xin12.counts) == "CNIH"] ="CNIH1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "SEPT4"] ="SEPTIN4"
rownames(Xin12.counts)[rownames(Xin12.counts) == "NEURL"] ="NEURL1"
rownames(Xin12.counts)[rownames(Xin12.counts) == "LARGE1"] ="LARGE"
rownames(Xin12.counts)[rownames(Xin12.counts) == "C7orf41"] ="MTURN"
Xin12.data <- CreateSeuratObject(counts = Xin12.counts, project = "Adult-X12")
original_dataset = rep("Adult_cadaveric_xin_sample_12", 2822)
age = rep(56,2822)
BMI = rep(21.2,2822)
gender = rep("male", 2822)
Xin12.data$original_dataset = original_dataset
Xin12.data$age = age
Xin12.data$BMI = BMI
Xin12.data$gender = gender

Xin.full <- merge(Xin1.data, y = c(Xin2.data,Xin3.data,Xin4.data,Xin5.data,Xin6.data,Xin7.data,Xin8.data,Xin9.data,Xin10.data,Xin12.data), add.cell.ids = c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X12"))

## Quality Control 
Xin.full[["percent.mt"]] = PercentageFeatureSet(Xin.full, pattern = "^MT-")
VlnPlot(Xin.full, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                            axis.title.y = element_text(size=0),
                                                                                            axis.text.x = element_text(size=0),
                                                                                            axis.text.y = element_text(size=20),
                                                                                            title = element_text(size=12))
Xin.QC <- subset(Xin.full, subset = nFeature_RNA > 1000 & nFeature_RNA < 3000 & percent.mt < 10 & nCount_RNA < 25000)

## PCA & Clustering
Xin.QC <- NormalizeData(Xin.QC)
Xin.QC <- FindVariableFeatures(Xin.QC, selection.method = "vst", nfeatures = 2000)
Xin.QC <- ScaleData(Xin.QC)
Xin.QC <- RunPCA(Xin.QC, features = VariableFeatures(object = Xin.QC))
Xin.QC <- FindNeighbors(Xin.QC, dims = 1:20)
Xin.QC <- FindClusters(Xin.QC, resolution = 0.7)
Xin.QC <- RunUMAP(Xin.QC, dims = 1:20)
DimPlot(Xin.QC, reduction = "umap", label = T)
FeaturePlot(Xin.QC, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
Xin.QC = RenameIdents(Xin.QC, '0' = "Alpha", '1' = "Beta",
                      '2' = "Alpha", '3' = "Beta", '4' = "CHGA-", 
                      '5' = "Delta", '6'="CHGA-", '7'="Alpha",
                      '8'= "CHGA-",'9'="PP", '10'="Beta", '11'="Beta",
                      '12'= "CHGA-",'13' = "CHGA-", '14'="Beta", '15' = "CHGA-",
                      '16' = "Prolif", '17' = "CHGA-")
Xin.QC$first_identity = Xin.QC@active.ident

## UMAP
DimPlot(Xin.QC, reduction = "umap", label = F, cols = c("#377EB8", "#E41A1C", "grey", "#984EA3", "#1B9E77", "#7570B3"))+NoLegend()

## Subsetting CHGA+ Cells
xin_endocrine = subset(x = Xin.QC, identity = c('Alpha', 'Beta', 'Delta', 'PP', 'Prolif'))


# Downloading Fetal human islet raw data from Cao et al. 2019 (DOI: 10.1126/science.aba7721)
# Downloaded from https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/

## Download Raw Data, Editing Gene Names, Creating Seurat Object, & Assigning Metadata
setwd("")
counts <- readRDS("")
cell.meta <- readRDS("")
gene.meta <- readRDS("")
gene.meta$gene_short_name <- as.character(gene.meta$gene_short_name)
gene.meta[gene.meta == "ATP5I"] <- "ATP5ME"
gene.meta[gene.meta == "ATP5G2"] <- "ATP5MC2"
gene.meta[gene.meta == "ATP5E"] <- "ATP5F1E"
gene.meta[gene.meta == "ATP5L"] <- "ATP5MG"
gene.meta[gene.meta == "TCEB2"] <- "ELOB"
gene.meta[gene.meta == "SHFM1"] <- "SEM1"
gene.meta[gene.meta == "USMG5"] <- "ATP5MD"
gene.meta[gene.meta == "NGFRAP1"] <- "BEX3"
gene.meta[gene.meta == "FAM213A"] <- "PRXL2A"
gene.meta[gene.meta == "GNB2L1"] <- "RACK1"
gene.meta[gene.meta == "MRP63"] <- "MRPL57"
gene.meta[gene.meta == "SF3B14"] <- "SF3B6"
gene.meta[gene.meta == "TCEB3"] <- "ELOA"
gene.meta[gene.meta == "C8orf59"] <- "RBIS"
gene.meta[gene.meta == "C19orf24"] <- "FAM174C"
gene.meta[gene.meta == "C20orf24"] <- "RAB5IF"
gene.meta[gene.meta == "FAM173A"] <- "ANTKMT"
gene.meta[gene.meta == "C6orf1"] <- "SMIM29"
gene.meta[gene.meta == "SEPT2"] <- "SEPTIN2"
gene.meta[gene.meta == "FAM57B"] <- "TLCD3B"
gene.meta[gene.meta == "PVRL2"] <- "NECTIN2"
gene.meta[gene.meta == "FAM96A"] <- "CIAO2A"
gene.meta[gene.meta == "KIAA1715"] <- "LNPK"
gene.meta[gene.meta == "RTFDC1"] <- "RTF2"
gene.meta[gene.meta == "FAM134A"] <- "RETREG2"
gene.meta[gene.meta == "STRA13"] <- "CENPX"
gene.meta[gene.meta == "C9orf3"] <- "AOPEP"
gene.meta[gene.meta == "PVRL3"] <- "NECTIN3"
gene.meta[gene.meta == "WRB"] <- "GET1"
gene.meta[gene.meta == "ASNA1"] <- "GET3"
gene.meta[gene.meta == "FAM60A"] <- "SINHCAF"
gene.meta[gene.meta == "FAM63B"] <- "MINDY2"
gene.meta[gene.meta == "COL4A3BP"] <- "CERT1"
gene.meta[gene.meta == "WHSC1"] <- "NSD2"
gene.meta[gene.meta == "APOPT1"] <- "COA8"
gene.meta[gene.meta == "TOMM70A"] <- "TOMM70"
gene.meta[gene.meta == "PVRL1"] <- "NECTIN1"
gene.meta[gene.meta == "TMEM261"] <- "DMAC1"
gene.meta[gene.meta == "FAM195A"] <- "MCRIP2"
gene.meta[gene.meta == "PPP2R4"] <- "PTPA"
gene.meta[gene.meta == "UQR11.1"] <- "UQCR11"
gene.meta[gene.meta == "C11orf73"] <- "HIKESHI"
gene.meta[gene.meta == "C19orf60"] <- "REX1BD"
gene.meta[gene.meta == "ZCCHC11"] <- "TUT4"
gene.meta[gene.meta == "LHFP"] <- "LHFPL6"
gene.meta[gene.meta == "MESDC2"] <- "MESD"
gene.meta[gene.meta == "FAM127B"] <- "RTL8A"
gene.meta[gene.meta == "C7orf55"] <- "FMC1"
gene.meta[gene.meta == "LINC00493"] <- "SMIM26"
gene.meta[gene.meta == "WBSCR22"] <- "BUD23"
gene.meta[gene.meta == "FAM105A"] <- "OTULINL"
gene.meta[gene.meta == "FAM127A"] <- "RTL8C"
gene.meta[gene.meta == "C14orf1"] <- "ERG28"
gene.meta[gene.meta == "SEPT2"] <- "SEPTIN2"
gene.meta[gene.meta == "VIMP"] <- "SELENOS"
gene.meta[gene.meta == "C16orf13"] <- "METTL26"
gene.meta[gene.meta == "SELM"] <- "SELENOM"
gene.meta[gene.meta == "FAM195B"] <- "MCRIP1"
gene.meta[gene.meta == "MLLT4"] <- "AFDN"
gene.meta[gene.meta == "WBP5"] <- "TCEAL9"
gene.meta[gene.meta == "SEPT11"] <- "SEPTIN11"
gene.meta[gene.meta == "C19orf70"] <- "MICOS13"
gene.meta[gene.meta == "WHSC1L1"] <- "NSD3"
gene.meta[gene.meta == "ATP5H"] <- "ATP5PD"
gene.meta[gene.meta == "TCEB1"] <- "ELOC"
gene.meta[gene.meta == "ATP5C1"] <- "ATP5F1C"
gene.meta[gene.meta == "SEP15"] <- "SELENOF"
gene.meta[gene.meta == "FAM96B"] <- "CIA02B"
gene.meta[gene.meta == "C17orf89"] <- "NDUFAF8"
gene.meta[gene.meta == "APOA1BP"] <- "NAXE"
gene.meta[gene.meta == "ATP5F1"] <- "ATP5PB"
gene.meta[gene.meta == "FAM46A"] <- "TENT5A"
gene.meta[gene.meta == "C14orf166"] <- "RTRAF"
gene.meta[gene.meta == "SEPT7"] <- "SEPTIN7"
gene.meta[gene.meta == "MINOS1"] <- "MICOS10"
gene.meta[gene.meta == "C11orf31"] <- "SELENOH"
gene.meta[gene.meta == "GLTSCR2"] <- "NOP53"
gene.meta[gene.meta == "LINC01420"] <- "NBDY"
gene.meta[gene.meta == "FAM159B"] <- "SHISAL2B"
gene.meta[gene.meta == "C7orf73"] <- "STMP1"
gene.meta[gene.meta == "C19orf43"] <- "TRIR"
gene.meta[gene.meta == "SELT"] <- "SELENOT"
gene.meta[gene.meta == "AES"] <- "TLE5"
gene.meta[gene.meta == "ATP5G1"] <- "ATP5MC1"
gene.meta[gene.meta == "ATP5B"] <- "ATP5F1B"
gene.meta[gene.meta == "SEPW1"] <- "SELENOW"
gene.meta[gene.meta == "ATP5D"] <- "ATP5F1D"
gene.meta[gene.meta == "ATP5O"] <- "ATP5PO"
gene.meta[gene.meta == "SELK"] <- "SELENOK"
gene.meta[gene.meta == "HN1"] <- "JPT1"
gene.meta[gene.meta == "ATPIF1"] <- "ATP5IF1"
gene.meta[gene.meta == "ATP5J"] <- "ATP5PF"
gene.meta[gene.meta == "ATP5A1"] <- "ATP5F1A"
gene.meta[gene.meta == "ATP5G3"] <- "ATP5MC3"
gene.meta[gene.meta == "ATP5J2"] <- "ATP5MF"
gene.meta[gene.meta == "MYEOV2"] <- "COPS9"
gene.meta[gene.meta == "C14orf2"] <- "ATP5MPL"
gene.meta[gene.meta == "UQCR11.1"] <- "UQCR11"
gene.meta[gene.meta == "C17orf76"] <- "SNHG29"
gene.meta[gene.meta == "SUZ12P"] <- "SUZ12P1"
gene.meta[gene.meta == "C19orf77"] <- "SMIM24"
gene.meta[gene.meta == "SEPT4"] <- "SEPTIN4"
gene.meta[gene.meta == "NEURL"] <- "NEURL1"
gene.meta[gene.meta == "LARGE1"] <- "LARGE"
gene.meta[gene.meta == "C7orf41"] <- "MTURN"
gene.meta$gene_short_name <- as.factor(gene.meta$gene_short_name)
rownames(counts) <- gene.meta$gene_short_name
rownames(cell.meta) <- cell.meta$sample
fetal.data <- CreateSeuratObject(counts = counts, project = "Fetal-C1", meta.data = cell.meta, min.cells = 3, min.features = 200)
fetal.data$sample = NULL
fetal.data$Exon_reads = NULL
fetal.data$Intron_reads = NULL
fetal.data$All_reads = NULL
fetal.data$RT_group = NULL
fetal.data$Organ = NULL
fetal.data$Batch = NULL
fetal.data$Main_cluster_umap_1 = NULL
fetal.data$Main_cluster_umap_2 = NULL
fetal.data$Main_cluster_name = NULL
fetal.data$Assay = NULL
fetal.data$Organ_cell_lineage = NULL
fetal.data$Experiment_batch = NULL
ids <- SplitObject(fetal.data, split.by = "Fetus_id")

fetus_1 = subset(fetal.data, subset = Fetus_id == "H27876")
original_dataset = rep("Fetal_cadaveric_cao_sample_1", 17984)
age = rep("Fetal",17984)
BMI = rep("Fetal",17984)
gender = rep(NA, 17984)
fetus_1$original_dataset = original_dataset
fetus_1$age = age
fetus_1$BMI = BMI
fetus_1$gender = gender

fetus_2 = subset(fetal.data, subset = Fetus_id == "H27870")
original_dataset = rep("Fetal_cadaveric_cao_sample_2", 30826)
age = rep("Fetal",30826)
BMI = rep("Fetal",30826)
gender = rep(NA, 30826)
fetus_2$original_dataset = original_dataset
fetus_2$age = age
fetus_2$BMI = BMI
fetus_2$gender = gender

fetus_3 = subset(fetal.data, subset = Fetus_id == "H27948")
original_dataset = rep("Fetal_cadaveric_cao_sample_3", 2828)
age = rep("Fetal",2828)
BMI = rep("Fetal",2828)
gender = rep(NA, 2828)
fetus_3$original_dataset = original_dataset
fetus_3$age = age
fetus_3$BMI = BMI
fetus_3$gender = gender

fetus_4 = subset(fetal.data, subset = Fetus_id == "H27464")
original_dataset = rep("Fetal_cadaveric_cao_sample_4", 775)
age = rep("Fetal",775)
BMI = rep("Fetal",775)
gender = rep(NA, 775)
fetus_4$original_dataset = original_dataset
fetus_4$age = age
fetus_4$BMI = BMI
fetus_4$gender = gender

fetus_5 = subset(fetal.data, subset = Fetus_id == "H27471")
original_dataset = rep("Fetal_cadaveric_cao_sample_5", 1648)
age = rep("Fetal",1648)
BMI = rep("Fetal",1648)
gender = rep(NA, 1648)
fetus_5$original_dataset = original_dataset
fetus_5$age = age
fetus_5$BMI = BMI
fetus_5$gender = gender

fetus_6 = subset(fetal.data, subset = Fetus_id == "H26547")
original_dataset = rep("Fetal_cadaveric_cao_sample_6", 3071)
age = rep("Fetal",3071)
BMI = rep("Fetal",3071)
gender = rep(NA, 3071)
fetus_6$original_dataset = original_dataset
fetus_6$age = age
fetus_6$BMI = BMI
fetus_6$gender = gender

fetus_7 = subset(fetal.data, subset = Fetus_id == "H26350")
original_dataset = rep("Fetal_cadaveric_cao_sample_7", 550)
age = rep("Fetal",550)
BMI = rep("Fetal",550)
gender = rep(NA, 550)
fetus_7$original_dataset = original_dataset
fetus_7$age = age
fetus_7$BMI = BMI
fetus_7$gender = gender

fetus_8 = subset(fetal.data, subset = Fetus_id == "H27058")
original_dataset = rep("Fetal_cadaveric_cao_sample_8", 2108)
age = rep("Fetal",2108)
BMI = rep("Fetal",2108)
gender = rep(NA, 2108)
fetus_8$original_dataset = original_dataset
fetus_8$age = age
fetus_8$BMI = BMI
fetus_8$gender = gender

fetus_9 = subset(fetal.data, subset = Fetus_id == "H27552")
original_dataset = rep("Fetal_cadaveric_cao_sample_9", 703)
age = rep("Fetal",703)
BMI = rep("Fetal",703)
gender = rep(NA, 703)
fetus_9$original_dataset = original_dataset
fetus_9$age = age
fetus_9$BMI = BMI
fetus_9$gender = gender

fetal.data= merge(fetus_1, c(fetus_2,fetus_3,fetus_4,fetus_5,fetus_6,fetus_7,fetus_8,fetus_9))

## Quality Control 
fetal.data[["percent.mt"]] = PercentageFeatureSet(fetal.data, pattern = "^MT-")
VlnPlot(fetal.data, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)& theme(axis.title.x = element_text(size=0),
                                                                                              axis.title.y = element_text(size=0),
                                                                                              axis.text.x = element_text(size=0),
                                                                                              axis.text.y = element_text(size=20),
                                                                                              title = element_text(size=12))
fetal.data <- subset(fetal.data, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & nCount_RNA < 5000)

## PCA & Clustering
fetal.data <- NormalizeData(fetal.data)
fetal.data <- FindVariableFeatures(fetal.data, selection.method = "vst", nfeatures = 2000)
fetal.data <- ScaleData(fetal.data)
fetal.data <- RunPCA(fetal.data, npcs = 30, verbose = FALSE)
fetal.data <- RunUMAP(fetal.data, reduction = "pca", dims = 1:30)
fetal.data <- FindNeighbors(fetal.data, reduction = "pca", dims = 1:30)
fetal.data <- FindClusters(fetal.data, resolution = 2.5)
DimPlot(fetal.data, reduction = "umap", label = T)
FeaturePlot(fetal.data, features = "CHGA", cols = brewer.pal(9, 'OrRd'))

## Assigning Cell Identity
fetal.data = RenameIdents(fetal.data, '0' = "CHGA-", '1' = "Poly",
                          '2' = "CHGA-", '3' = "CHGA-", '4' = "CHGA-", 
                          '5' = "CHGA-", '6'="CHGA-", '7'="CHGA-",
                          '8'="Beta",'9'="CHGA-", '10'="CHGA-", '11'="CHGA-",
                          '12'="CHGA-",'13' = "CHGA-", '14'="Alpha", '15'="CHGA-",
                          '16'="CHGA-",'17'= "CHGA-", '18' = "CHGA-", '19' = "CHGA-", '20' = "CHGA-"
                          , '21' = "Beta", '22' = "Delta", '23' = "CHGA-", '24' = "CHGA-"
                          , '25' = "CHGA-", '26' = "CHGA-", '27' = "CHGA-", '28' = "CHGA-"
                          , '29' = "CHGA-", '30' = "CHGA-", '31' = "CHGA-", '32' = "CHGA-", '33' = "Delta"
                          , '34' = "NeuroEndo", '35' = "CHGA-", '36' = "Beta", '37' = "CHGA-", '38' = "CHGA-", '39' = "CHGA-"
                          , '40' = "CHGA-", '41' = "CHGA-", '42' = "CHGA-", '43' = "CHGA-", '44' = "CHGA-", '44' = "CHGA-", '45' = "CHGA-"
                          , '46' = "CHGA-", '47' = "CHGA-", '48' = "CHGA-", '49' = "CHGA-", '50' = "CHGA-", '51' = "Beta", '52' = "CHGA-"
                          , '53' = "CHGA-", '54' = "Beta", '55' = "CHGA-", '56' = "CHGA-", '57' = "CHGA-", '58' = "NeuroEndo")
fetal.data$CHGA_identity = fetal.data@active.ident

## UMAP
DimPlot(fetal.data, reduction = "umap", label = F, cols = c("grey", "#FF7F00", "3E41A1C", "#377EB8", "#984EA3", "#A6761D"))+NoLegend()

## Subsetting CHGA+ Cells
fetal_endocrine = subset(x = fetal.data, identity = c('Poly', 'Beta', 'Alpha', 'Delta', 'NeuroEndo'))


# Merging All CHGA+ Populations
all_CHGA = merge(augs_endocrine, c(veres_endocrine, weng_endocrine, balboa_endocrine,
                                   balboa_vivo_1_endocrine, balboa_vivo_6_endocrine, 
                                   augs_vivo_6_endocrine, baron_endocrine, fang_endocrine,
                                   xin_endocrine, fetal_endocrine))


# Addition of Sequencing Metadata 
added_column = c()
length = ncol(all_CHGA)
count = 1
while (count <= length) {
  if(all_CHGA$orig.ident[count] %in% c("Augsornworawat")){
    value = "Aug SC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("veres")){
    value = "Veres SC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("weng")){
    value = "Weng SC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Balboa")){
    value = "Balboa SC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Balboa vivo 1 month")){
    value = "Balboa Vivo 1 month"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Balboa vivo 6 month")){
    value = "Balboa Vivo 6 month"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Augs vivo")){
    value = "Aug Vivo 6 month"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Adult-F1","Adult-F2","Adult-F3","Adult-F4","Adult-F6")){
    value = "Fang"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Adult-B1","Adult-B2","Adult-B3","Adult-B4")){
    value = "Baron"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Adult-X1","Adult-X2","Adult-X3","Adult-X4","Adult-X5","Adult-X6","Adult-X7","Adult-X8","Adult-X9","Adult-X10","Adult-X12")){
    value = "Xin"
  }
  else{
    value = "Cao"
  }
  added_column = append(added_column, value)
  count = count + 1
}
View(added_column)
all_CHGA$sequencing = added_column

## Addition of Sequencing2 Column 
added_column = c()
length = ncol(all_CHGA)
count = 1
while (count <= length) {
  if(all_CHGA$orig.ident[count] %in% c("Augsornworawat", "Augs vivo")){
    value = "Aug SC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("veres")){
    value = "Veres SC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("weng")){
    value = "Weng SC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Balboa", "Balboa vivo 1 month", "Balboa vivo 6 month")){
    value = "Balboa SC"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Fetal-C1")){
    value = "Fetal"
  }
  else{
    value = "Adult"
  }
  added_column = append(added_column, value)
  count = count + 1
}

View(added_column)
all_CHGA$sequencing2 = added_column

## Addition of Maturation Metadata
added_column = c()
length = ncol(all_CHGA)
count = 1
while (count <= length) {
  if(all_CHGA$orig.ident[count] %in% c("Augsornworawat", "veres", "weng", "Balboa")){
    value = "SC-Vitro"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Balboa vivo 1 month", "Balboa vivo 6 month", 
                                            "Augs vivo")){
    value = "SC-Vivo"
  }
  else if(all_CHGA$orig.ident[count] %in% c("Fetal-C1")){
    value = "Fetal"
  }
  else{
    value = "Adult"
  }
  added_column = append(added_column, value)
  count = count + 1
}
View(added_column)
all_CHGA$maturation = added_column


# Removing all mitochondrial genes (MT) from all data for consistency 
DefaultAssay(all_CHGA) = "RNA"
all_genes = rownames(all_CHGA)
anti_mt_genes = all_genes[grep(x = all_genes, pattern ="^MT-", invert = T)]
mt_genes = all_genes[grepl(x = all_genes, pattern = "MT-")]
all_CHGA = subset(all_CHGA, features = anti_mt_genes)

# Setting Up Object for Integration
split.list <- SplitObject(all_CHGA, split.by = "sequencing")
split.list <- lapply(X = split.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = split.list)

# Performing Integration
endocrine.anchors <- FindIntegrationAnchors(object.list = split.list, anchor.features = features)
CHGA_integrated <- IntegrateData(anchorset = endocrine.anchors)
DefaultAssay(CHGA_integrated) <- "integrated"
CHGA_integrated <- ScaleData(CHGA_integrated, vars.to.regress = c('percent.mt', 'nCount_RNA'))
CHGA_integrated <- RunPCA(CHGA_integrated, npcs = 30, verbose = FALSE)
ElbowPlot(CHGA_integrated)
CHGA_integrated <- RunUMAP(CHGA_integrated, reduction = "pca", dims = 1:15)
CHGA_integrated <- FindNeighbors(CHGA_integrated, reduction = "pca", dims = 1:15)
CHGA_integrated <- FindClusters(CHGA_integrated, resolution = 2)

# Cell Typing
CHGA_integrated = RenameIdents(CHGA_integrated, '0' = "Beta", '1' = "Alpha",
                               '2' = "Alpha", '3' = "Alpha", '4' = "Alpha", 
                               '5' = "EC", '6'="Beta", '7'="Delta",
                               '8'="Alpha",'9'="Beta", '10'="Beta", '11'="Alpha",
                               '12'="Beta",'13' = "Beta", '14'="Beta", '15'="EC",
                               '16'="Alpha",'17'= "Beta", '18' = "Acinar", '19' = "Poly", '20' = "PP"
                               , '21' = "Prolif", '22' = "Acinar", '23' = "NeuroEndo", '24' = "Poly"
                               , '25' = "EC", '26' = "Alpha", '27' = "Beta", '28' = "Acinar"
                               , '29' = "Beta", '30' = "Prog", '31' = "Beta", '32' = "Acinar", '33' = "Delta"
                               , '34' = "Beta", '35' = "Epsilon")
CHGA_integrated$CHGA_identity = CHGA_integrated@active.ident

# Removing CHGA+ Acinar Cells to generate final endocrine integrated dataset
endocrine_integrated = subset(CHGA_integrated, subset = CHGA_identity %in% c("Beta", "Alpha", "Delta", "PP", "EC", "Epsilon", "Prolif", "Prog", 'Poly', 'NeuroEndo'))
saveRDS(endocrine_integrated, ".RDS")

