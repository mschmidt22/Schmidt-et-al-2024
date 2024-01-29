# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig S4

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

## Dotplot of genes involved in INS secretion (Fig S4A)
## Gene set from JAX (http://www.informatics.jax.org/go/term/GO:0032024#myDataTable=results%3D100%26startIndex%3D0%26sort%3Dterm%26dir%3Dasc)
pos_reg_INS_secr_GO = c('AACS', 'ABAT', 'ABCC8', 'ANXA7', 'BAD', 'CACNA1C', 'CACNA1D', 'CAMK2N1', 'CASK', 'CASR',
                        'DYNLL1', 'GCK', 'GIPR', 'GLP1R', 'GLUD1', 'GLUL', 'GNA11', 'GNAQ', 'GNAS', 'HIF1A',
                        'IRS2', 'ISL1', 'KIF5B', 'MPC2', 'NKX6-1', 'NNAT', 'OGA', 'OSBP', 'OXCT1', 'PDX1', 'PFKFB2',
                        'PFKM', 'PPP3CB', 'PRKACA', 'PRKAR1A', 'PSMD9', 'PTBP1', 'RAC1', 'RBP4', 'RFX6', 'RPH3AL',
                        'SERP1', 'SLC30A8', 'SOX4', 'SRI', 'TARDBP', 'TM7SF3', 'UCN3')
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
DotPlot(all_beta, features = pos_reg_INS_secr_GO, cols = c("lightgrey", "red"), dot.scale = 20)+RotatedAxis()+theme(axis.text.x = element_text(size = 20))

## Dotplot of genes involved in Beta Identity (Fig S4A)
## Loading Gene Set from van Gurp et al. 2020 (DOI: 10.1038/s41467-022-29588-8)
beta_identity_genes = c('ABCC8', 'ADCYAP1', 'ALCAM', 'ARL6IP5', 'ATP2A3', 'ATP6V1A', 'C1QL1', 'CADM1', 'CASR', 'CDKN1C',
                        'CNP', 'DHRS2', 'DHRS7', 'DLK1', 'DNAJB9', 'EDARADD', 'EIF4A2','ELAVL4', 'ENO1', 'ENTPD3', 'ERO1B',
                        'ETNK1', 'OTULINL', 'SHISAL2B', 'FXYD2', 'G3BP1', 'G6PC2', 'GAD2', 'GLIS3', 'GSN', 'HADH',
                        'HSP90AB1', 'IAPP', 'IGF1R', 'IGF2', 'INS', 'MAFA', 'MAP1B', 'MGAT4A', 'MXRA7', 'NKX6-1',
                        'NPTX2', 'NTPCR', 'PCSK1', 'PDX1', 'PEBP1', 'PFN2', 'PLCB4', 'PLCXD3', 'PLEKHA6', 'PRSS23',
                        'PRUNE2', 'PTEN', 'RAP1GAP2', 'RBP4', 'ROBO2', 'RPL3', 'RPS3', 'SCD', 'SCG3', 'SCGN', 'SERINC1',
                        'SLC30A8', 'SORL1', 'SRPRB', 'STX1A', 'SURF4', 'SUSD4', 'SYNE2', 'SYT13', 'TIMP2', 'TMOD1',
                        'TPM3', 'TSPAN1', 'TSPAN13', 'TSPYL1', 'UCHL1')
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
DotPlot(all_beta, features = beta_identity_genes, cols = c("lightgrey", "red"), dot.scale = 20)+RotatedAxis()+theme(axis.text.x = element_text(size = 20))

