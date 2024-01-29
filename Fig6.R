# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig 6

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

## Percent of Cells Expressing Shared Pancreas-Neuron TFs (Fig 6A)
neuron_TFs = c('ASCL1', 'DNMT3A', 'ETV1', 'FOXA1', 'FOXA2', 'HDAC1' ,'HDAC2', 'HIPK2', 'ISL1',
               'MEIS2', 'NEUROD1', 'NKX2-2', 'NKX6-1', 'ONECUT2', 'PAX6' , 'PBX1', 'PROX1', 
               'RBPJ', 'SMAD2', 'SOX11', 'SOX4')
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
SC_vitro_neural_pan_DE_TFs = FindMarkers(all_beta, ident.1 = "SC-Vitro", ident.2 = c("Adult", "Fetal", 'SC-Vivo'), features = neuron_TFs, logfc.threshold = 0, min.pct = 0)
SC_vivo_neural_pan_DE_TFs = FindMarkers(all_beta, ident.1 = "SC-Vivo", ident.2 = c("SC-Vitro", "Adult", 'Fetal'), features = neuron_TFs, logfc.threshold = 0, min.pct = 0)
adult_neural_pan_DE_TFs = FindMarkers(all_beta, ident.1 = "Adult", ident.2 = c("SC-Vitro", "Fetal", 'SC-Vivo'), features = neuron_TFs, logfc.threshold = 0, min.pct = 0)
fetal_neural_pan_DE_TFs = FindMarkers(all_beta, ident.1 = "Fetal", ident.2 = c("SC-Vitro", "Adult", 'SC-Vivo'), features = neuron_TFs, logfc.threshold = 0, min.pct = 0)

## Heatmap of shared neural-pancreatic target genes (Fig 6B)

### Gene sets identified from SCENIC analysis target genes & cross referenced with GOBP_NEURON_DEVELOPMENT
### https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_NEURON_DEVELOPMENT.html
pbx1_neural_targets = c('BASP1', 'CDH2', 'CELF4', 'CLSTN3', 'FZD3', 'GABRB3', 
                          'KIF1B', 'KIF5A', 'NF1', 'PSD3', 'RAB3B', 'RIMS2', 'SV2A', 'USP9X')
sox4_neural_targets = c('ACAP3', 'APOE', 'CADPS', 'CDH8', 'CXADR', 'DYNC1H1', 'NEFL', 
                          'NIPSNAP1', 'NSG1', 'PRKAR1A', 'PTPRN2', 'RPS6KA3', 'STRN4', 'YWHAG')
pax4_neural_targets = c('C1QL1', 'CNIH2', 'CTNNB1', 'GAP43', 'GLUL', 'HNRNPH2', 'KIF1A',
                          'MAL2', 'NEFM', 'PALLD', 'RAB3B', 'RTN1', 'SCAMP5', 'SYT11')
isl1_neural_targets = c('ATP2B4', 'BAALC', 'FYN', 'KLC1', 'MAP2', 'NFASC', 'PAK2', 'RAB5B',
                          'RUFY3', 'SEMA3A', 'STX3', 'SYNJ2', 'UTRN')
sox11_neural_targets = c('DBN1', 'DCC', 'DDC', 'ELOVL5', 'FARP1', 'FXYD6', 'GNAI2', 'GNAQ',
                           'IGF2BP1', 'INA', 'ITGB1', 'MAP1B', 'MAP6', 'STMN2')
smad1_neural_targets = c('ARL8A', 'ARMCX3', 'CRKL', 'CTNND2', 'DPYSL5', 'DST', 'DTNA', 
                           'INPP5F', 'KIF5C', 'KLHL24', 'SLC17A6', 'STMN4', 'TUBA1A', 'YWHAZ')
nkx22_neural_targets = c('ADNP', 'AKAP9', 'AP2M1', 'AP2S1', 'AUTS2', 'FBXW11', 'KCTD16',
                           'MOB4', 'PRKCI', 'PRR7', 'RAB3D', 'RPS17', 'SRSF10', 'STAU1')
dnmt3a_neural_targets = c('ATP2A2', 'CAPRIN1', 'CLCN3', 'CPEB4', 'DCTN1', 'EIF3B', 'EIF4A3',
                            'ENAH', 'GNB1', 'GRIPAP1', 'IGF1R', 'MPST', 'RAB11A', 'TTLL7')


### Scaling Data
DefaultAssay(all_beta) = 'RNA'
all_genes = rownames(all_beta)
all_beta = ScaleData(all_beta, features = all_genes)
all_beta_shuffle = RandomSubsetData(all_beta, .9999)
DefaultAssay(all_beta_shuffle) = "RNA"
Idents(all_beta_shuffle) = "maturation"
my_levels <- c('SC-Vitro', 'SC-Vivo', 'Adult', 'Fetal')
all_beta_shuffle@active.ident <- factor(x = all_beta_shuffle@active.ident, levels = my_levels)

### Heatmaps
DoHeatmap(shuffle, features = pbx1_neural_targets, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
DoHeatmap(shuffle, features = sox4_neural_targets, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
DoHeatmap(shuffle, features = pax4_neural_targets, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
DoHeatmap(shuffle, features = isl1_neural_targets, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
DoHeatmap(shuffle, features = sox11_neural_targets, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
DoHeatmap(shuffle, features = smad1_neural_targets, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
DoHeatmap(shuffle, features = nkx22_neural_targets, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
DoHeatmap(shuffle, features = dnmt3a_neural_targets, slot = 'scale.data')+ scale_fill_gradientn(colors = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))

## Figure 6c was generated by running the same SCENIC Analysis as in Fig4
## Transcription factors the highest AUC were then cross referenced if any of their targets were the same genes that were previvously identified in Fig5D
