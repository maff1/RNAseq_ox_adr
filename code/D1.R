#!/usr/bin/env Rscript
################################################################################
## RNAseq workflow - PART D:: visual exploration ###############################
################################################################################
# packages <- c("genefilter", "biomaRt", "pheatmap", "png",
#               "dplyr", "RColorBrewer", "ggplot2", "EnhancedVolcano")
# ipak <- function(pkg) {
#     new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#     if (length(new.pkg))
#         install.packages(new.pkg, dependencies = TRUE,
#                          repos = "http://cran.r-project.org/")
#     sapply(pkg, require, character.only = TRUE)
# }
# ipak(packages)
# BiocManager::install("EnhancedVolcano")
################################################################################
library(ggplot2)
library(dplyr)
library(genefilter)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(png)
library(EnhancedVolcano)
library(ggrepel)
################################################################################
## HEATMAP #####################################################################
topVarGenes <- head(order(-rowVars(assay(vsd))),20)

geneNames <- ensembldb::mapIds(x = EDB,
                               keys = rownames(assay(vsd))[topVarGenes],
                               column = "SYMBOL",
                               keytype = "GENEID",
                               multiVals = "first")
ann_colors = list(
    Time = c("white", "firebrick"),
    CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
    GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)
df <- as.data.frame(colData(dds.un.shr)[,c("Sample","Genotype", "Condition")])

color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                            "RdYlBu")))(100)
heatm.vsd <- pheatmap(assay(vsd)[topVarGenes,], labels_row=geneNames,
                      cluster_rows=T, show_rownames=T, show_colnames=F,
                      cluster_cols=T, annotation_col=df,
                      #annotation_colors = brewer.pal(n = 5, name ="Set3"),
                      cutree_rows = 2, cutree_cols = 2,
                      clustering_distance_rows="euclidean",
                      clustering_distance_cols="euclidean",
                      clustering_method="complete", border_color=FALSE,
                      scale="row", cex=1,
                      color=color,
                      main="top20 var. genes")
ggsave("heatmapMod_exp2.jpg", plot = heatm.vsd, device = NULL,
       path = "./results",
       scale = 1, width = 20, height = 20, units = "cm",
       dpi = 600)
