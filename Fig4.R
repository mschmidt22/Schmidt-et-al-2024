# Schmidt et al. 2024 (DOI: 10.1186/s12864-024-10013-x)
# Fig 4

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
sc_adult_beta = subset(all_beta, maturation %in% c("SC-Vitro", "Adult"))
adult_fetal_beta = subset(all_beta, maturation %in% c("Adult", "Fetal"))
sc_adult_fetal_beta = subset(all_beta, maturation %in% c("SC-Vitro", "Adult", "Fetal"))

## Downloading Transcription Factor List
## Adopted from Lambert et al. 2018 (DOI: 10.1016/j.cell.2018.01.029)
TF = read.csv("tf database.txt", header = F, sep = ",", dec = ".")
TF <- TF %>% separate(V1,c("TF","Target", "Mechanism", "Number"),remove=T)
TF_list = TF$TF
TF_list = unique(TF_list)

## Downloading Pancreatic Development Gene List
## (http://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_PANCREAS_DEVELOPMENT)
pancreas_dev_genes = read.delim('GOBP_PANCREAS_DEVELOPMENT.v2022.1.Hs.txt')
pancreas_dev_genes = pancreas_dev_genes$GOBP_PANCREAS_DEVELOPMENT[2:83]

## SC vs Adult Beta DE TFs Volcano Plot (Fig 4A)
DefaultAssay(sc_adult_beta) = "RNA"
Idents(sc_adult_beta) = "maturation"
SC_vitro_adult_beta_DEGs = FindAllMarkers(sc_adult_beta, logfc.threshold = 0)
DEG_genes = rownames(SC_vitro_adult_beta_DEGs)
DEG_and_TFs = append(DEG_genes, TF_list)
beta_DEG_TFs = DEG_and_TFs[duplicated(DEG_and_TFs)]
beta_DEG_TFs = append(beta_DEG_TFs, c("NKX2-2", 'NKX6-1', 'INSM1'))
DEG_TFs = FindMarkers(sc_adult_beta, ident.1 = "Adult", ident.2 = "SC-Vitro", features = beta_DEG_TFs, min.pct = 0, min.diff.pct = 0, logfc.threshold = 0)
DEG_TFs = subset(DEG_TFs, p_val_adj<0.05)
EnhancedVolcano(DEG_TFs, lab = rownames(DEG_TFs),
                x = 'avg_log2FC', y = 'p_val_adj',xlim = c(-2.5, 2.55),title = 'Adult vs SC',
                pCutoff = NA,FCcutoff = NA, pointSize = 3,labSize = 0,colAlpha = 1.0,
                col=c('grey', 'black','grey'),gridlines.minor = FALSE,
                gridlines.major = FALSE, legendPosition = 'top',drawConnectors = F,widthConnectors = 0,
                colConnectors = 'black', cutoffLineCol = "red")

## Adult vs Fetal Beta DE TFs Volcano Plot (Fig 4B)
DefaultAssay(adult_fetal_beta) = "RNA"
Idents(adult_fetal_beta) = "maturation"
adult_fetal_beta_DEGs = FindAllMarkers(adult_fetal_beta, logfc.threshold = 0)
DEG_genes_2 = rownames(adult_fetal_beta_DEGs)
DEG_and_TFs_2 = append(DEG_genes_2, TF_list)
beta_DEG_TFs_2 = DEG_and_TFs_2[duplicated(DEG_and_TFs_2)]
beta_DEG_TFs_2 = append(beta_DEG_TFs_2, c("NKX2-2", 'NKX6-1', 'INSM1'))
DEG_TFs_2 = FindMarkers(adult_fetal_beta, ident.1 = "Adult", ident.2 = "Fetal", features = beta_DEG_TFs_2, min.pct = 0, min.diff.pct = 0, logfc.threshold = 0)
DEG_TFs_2 = subset(DEG_TFs_2, p_val_adj<0.05)
EnhancedVolcano(DEG_TFs_2, lab = rownames(DEG_TFs_2),
                x = 'avg_log2FC', y = 'p_val_adj',xlim = c(-2.5, 2.5),title = 'Adult vs Fetal',
                pCutoff = NA,FCcutoff = NA, pointSize = 3,labSize = 0,colAlpha = 1.0,
                col=c('grey', 'black','grey'),gridlines.minor = FALSE,
                gridlines.major = FALSE, legendPosition = 'top',drawConnectors = F,widthConnectors = 0,
                colConnectors = 'black', cutoffLineCol = "red")

## Percent of Cells Expressing Pancreatic Development Relevant TFs (Fig 4c)
pancreas_and_TFs = append(pancreas_dev_genes, TF_list)
pancreas_TFs = pancreas_and_TFs[duplicated(pancreas_and_TFs)]
DefaultAssay(all_beta) = "RNA"
Idents(all_beta) = "maturation"
SC_vitro_pancreatic_TF_expression = FindMarkers(all_beta, ident.1 = "SC-Vitro", ident.2 = c("Adult", "Fetal", 'SC-Vivo'), features = c(pancreas_TFs, "NKX6-1", "NKX2-2", "PBX1"), logfc.threshold = 0, min.pct = 0)
SC_vivo_pancreatic_TF_expression = FindMarkers(all_beta, ident.1 = "SC-Vivo", ident.2 = c("SC-Vitro", "Adult", 'Fetal'), features = c(pancreas_TFs, "NKX6-1", "NKX2-2", "PBX1"), logfc.threshold = 0, min.pct = 0)
adult_pancreatic_TF_expression = FindMarkers(all_beta, ident.1 = "Adult", ident.2 = c("SC-Vitro", "Fetal", 'SC-Vivo'), features = c(pancreas_TFs, "NKX6-1", "NKX2-2", "PBX1"), logfc.threshold = 0, min.pct = 0)
fetal_pancreatic_TF_expression = FindMarkers(all_beta, ident.1 = "Fetal", ident.2 = c("SC-Vitro", "Adult", 'SC-Vivo'), features = c(pancreas_TFs, "NKX6-1", "NKX2-2", "PBX1"), logfc.threshold = 0, min.pct = 0)

## Regulon Analysis (Fig 4D-G)

### Setting up regulon 
build_loom(file.name = "input.loom",
           dgem = sc_adult_fetal_beta@assays$RNA@counts,
           title = "SC, Adult, Fetal Beta Comparison",
           default.embedding = sc_adult_fetal_beta@reductions$umap@cell.embeddings,
           default.embedding.name = "UMAP")
loom <- open_loom("input.loom", mode = "r+")
add_col_attr(loom = loom, key = "orig.ident", value = sc_adult_fetal_beta@meta.data$orig.ident)
add_col_attr(loom = loom, key = "cell_label", value = sc_adult_fetal_beta@meta.data$cell_label)
add_col_attr(loom = loom, key = "age", value = sc_adult_fetal_beta@meta.data$age)
add_col_attr(loom = loom, key = "sequencing", value = sc_adult_fetal_beta@meta.data$sequencing)
add_col_attr(loom = loom, key = "maturation", value = sc_adult_fetal_beta@meta.data$maturation)
add_col_attr(loom = loom, key = "endocrine_identity", value = sc_adult_fetal_beta@meta.data$endocrine_identity)
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
sc_adult_fetal_beta[["AUC"]] <- CreateAssayObject(data = AUCmat)
DefaultAssay(sc_adult_fetal_beta) <- "AUC"
sc_adult_fetal_beta <- ScaleData(sc_adult_fetal_beta, assay = "AUC")

### Calculate & Plotting RSS S-plots (Fig 4D)
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
                     box.padding   = .35, 
                     point.padding = .5,
                     segment.color = 'black',
                     na.rm=TRUE,
                     segment.alpha=1,
                     alpha=1,
                     seed = 555)+
    theme_classic()
}
plotRSS_oneSet(rss, setName = "SC-Vitro", n =15)
plotRSS_oneSet <- function(rss, setName, n=15)
{
  library(ggplot2)
  library(ggrepel)
  
  rssThisType <- sort(rss[,setName], decreasing=TRUE)
  thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
  thisRss$regulon[(n+1):nrow(thisRss)] <- NA
  
  ggplot(thisRss, aes(x=rank, y=rss)) + 
    ggtitle(setName) + 
    geom_point(color='dodgerblue',size=3)+
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
plotRSS_oneSet(rss, setName = "Adult", n =15)
plotRSS_oneSet <- function(rss, setName, n=15)
{
  library(ggplot2)
  library(ggrepel)
  
  rssThisType <- sort(rss[,setName], decreasing=TRUE)
  thisRss <- data.frame(regulon=names(rssThisType), rank=seq_along(rssThisType), rss=rssThisType)
  thisRss$regulon[(n+1):nrow(thisRss)] <- NA
  
  ggplot(thisRss, aes(x=rank, y=rss)) + 
    ggtitle(setName) + 
    geom_point(color='green3',size=3)+
    geom_label_repel(size = 5, aes(label = regulon),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'black',
                     na.rm=TRUE,
                     segment.alpha=1,
                     alpha=1,
                     seed = 555)+
    theme_classic()
}
plotRSS_oneSet(rss, setName = "Fetal", n =15)

### RSS Heatmap (Fig 4E)
selectedResolution <- 'maturation'
cells_per_maturation <- split(rownames(maturation), maturation[,selectedResolution]) 
regulonActivity_by_maturation <- sapply(cells_per_maturation,
                                        function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_by_maturation_scaled <- t(scale(t(regulonActivity_by_maturation), center = T, scale=T))
write.csv(regulonActivity_by_maturation_scaled, 'regulon_scaled.csv')
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_by_maturation_scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6)))

### Regulon UMAP (Fig 4F)
DefaultAssay(sc_adult_fetal_beta) <- "AUC"
sc_adult_fetal_beta <- FindVariableFeatures(sc_adult_fetal_beta, selection.method = "vst", nfeatures = 2000)
sc_adult_fetal_beta <- ScaleData(sc_adult_fetal_beta, assay = "AUC", vars.to.regress = c('sequencing'))
sc_adult_fetal_beta <- RunPCA(sc_adult_fetal_beta, features = VariableFeatures(sc_adult_fetal_beta))
sc_adult_fetal_beta <- FindNeighbors(sc_adult_fetal_beta, dims = 1:20)
sc_adult_fetal_beta <- FindClusters(sc_adult_fetal_beta, resolution = 0.5)
sc_adult_fetal_beta <- RunUMAP(sc_adult_fetal_beta, dims = 1:20)
DimPlot(sc_adult_fetal_beta, reduction = "umap", label = F, group.by = 'maturation', cols = c('dodgerblue', 'green3', 'magenta2'))+NoLegend()

### Expression & Activity Featureplots (Fig 4G) 
DefaultAssay(sc_adult_fetal_beta) = "RNA"
FeaturePlot(sc_adult_fetal_beta, features = "SOX4")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(sc_adult_fetal_beta, features = "SOX11")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(sc_adult_fetal_beta, features = "PBX1")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
FeaturePlot(sc_adult_fetal_beta, features = "PAX4")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+NoLegend()
DefaultAssay(sc_adult_fetal_beta) = "AUC"
FeaturePlot(sc_adult_fetal_beta, features = "SOX4(+)")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+NoLegend()
FeaturePlot(sc_adult_fetal_beta, features = "SOX11(+)")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+NoLegend()
FeaturePlot(sc_adult_fetal_beta, features = "PBX1(+)")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+NoLegend()
FeaturePlot(sc_adult_fetal_beta, features = "PAX4(+)")+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))+NoLegend()

