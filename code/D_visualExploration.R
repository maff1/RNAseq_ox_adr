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
# mat <- assay(rld)[ topVarGenes, ]
# mat <- mat - rowMeans(mat)
# df <- as.data.frame(colData(rld)[,c("cell","dex")])
# pheatmap(mat, annotation_col=df)

# ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
# genemap <- getBM( attributes = c("ensembl_gene_id","hgnc_symbol"),
#                   filters = "ensembl_gene_id",
#                   values = rownames(assay(vsd)),
#                   mart = ensembl)

geneNames <- ensembldb::mapIds(x = EDB,
                      keys = rownames(assay(vsd))[topVarGenes],
                      column = "SYMBOL",
                      keytype = "GENEID",
                      multiVals = "first")


# geneDat <- rownames(assay(vsd)[topVarGenes,])
# idx <- match(geneDat, genemap$ensembl_gene_id)
# geneNames <- genemap$hgnc_symbol[idx]

# select <- order(rowMeans(counts(ddsSeq,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(ddsSeq)[,c("Phenotype","Treatment", "Condition")])
# heatm <- pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=TRUE,
#          cluster_cols=T, annotation_col=df)
# ggsave("heatmap.jpg", plot = heatm, device = NULL,
#path = "F:/data/RNAseq/P190850/analysis/redo",
#        scale = 0.8, width = 20, height = 20, units = "cm",
#        dpi = 600)

# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:100]

# rownames(assay(vsd)[topVarGenes,]) <- geneNames
# Specify colors
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

# Sample to sample distances ###################################################
# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)), method="euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Sample
colnames(sampleDistMatrix) <- vsd$Sample
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatm.ss <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         clustering_method="complete",
         cex=1,
         col=colors, border_color=FALSE,
         cutree_rows = 2, cutree_cols = 2,
         main="sample-to-sample distances",
         show_colnames=T)
################################################################################
## PCA #########################################################################
pcaData2 <- plotPCA(vsd, intgroup=c("Genotype", "Condition"), returnData=TRUE, ntop = nrow(vsd))
percentVar <- round(100 * attr(pcaData2, "percentVar"))
pcaData2$Sample <- df$Sample
mypca2 <- ggplot(pcaData2, aes(PC1, PC2, color=Genotype, shape= Condition,
                               #label=gsub("WTCHG_799608_","",pcaData2$name))) +
                               label=Sample)) +
    geom_point(size=3) +
    geom_label_repel(show.legend = FALSE) +
    xlim(-100, 100) + ylim(-100, 100) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_bw() + theme(axis.text.x = element_text(color="black", size=12),
                       axis.text.y = element_text(color="black", size=12),
                       axis.title.x = element_text(color="black", size=14, face="bold"),
                       axis.title.y = element_text(color = "black", size=14, face="bold"),
                       legend.text = element_text(colour="black", size = 12, face = "bold"),
                       legend.title = element_text(colour="blue", size = 12, face = "bold"),
                       plot.title = element_text(hjust = 0.5, color = "black", size=14, face="bold")) +
    ggtitle("PCA - scorings") +
geom_hline(yintercept=0, linetype="dashed",
            color = "black", size=0.5)+
geom_vline(xintercept=0, linetype="dashed",
            color = "black", size=0.5)
ggsave("pcaMod_EXP2.jpg", plot = mypca2,
       path = "./results",
       scale = 1, width = 20, height = 20, units = "cm",
       dpi = 600)
#############################################################################
mastyle <- theme_bw() +
    theme(axis.text.x = element_text(color="black", size=12),
                              axis.text.y = element_text(color="black", size=12),
                              axis.title.x = element_text(color="black", size=14, face="bold"),
                              axis.title.y = element_text(color = "black", size=14, face="bold"),
                              legend.text = element_text(colour="black", size = 12, face = "bold"),
                              legend.title = element_text(colour="blue", size = 12, face = "bold"),
                            plot.title = element_text(hjust = 0.5, color = "black", size=14, face="bold"))
#############################################################################
library(factoextra)
library(FactoMineR)

mat <- t(assay(vsd))
mat2 <- mat
colnames(mat2) <- ensembldb::mapIds(x = EDB,
                               keys = colnames(mat2),
                               column = "SYMBOL",
                               keytype = "GENEID",
                               multiVals = "first")

#table(duplicated(colnames(mat2)))
#column_duped <- apply(mat2, MARGIN = 2, duplicated)
myres <- mat2[, !duplicated(colnames(mat2))]
res.pca <- prcomp(mat, scale = F)
res.pca2 <- prcomp(myres, scale = F)
samples <- as.factor(df$Sample)
biplot <- fviz_pca_biplot(res.pca, select.var=list(contrib=10), repel = TRUE,
                addEllipses = F,
                col.ind = samples,
                pointsize=1, pointshape=21, fill.ind = samples,
                legend.title = "samples",
                title = "PCA - Biplot") + xlim(-100, 100) + ylim (-100, 100) + mastyle


ggsave("biplot_EXP2.jpg", plot = biplot,
       path = "./results",
       scale = 1, width = 20, height = 20, units = "cm",
       dpi = 600)

pca.load <- fviz_pca_var(res.pca, select.var=list(contrib=20), repel = TRUE,
             addEllipses = F, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             title = "PCA - loadings") + xlim(-4, 4) + ylim (-4, 4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + mastyle

pca.load2 <- fviz_pca_var(res.pca2, select.var=list(contrib=20), repel = TRUE,
                         addEllipses = F, col.var = "contrib",
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         title = "PCA - loadings") + xlim(-4, 4) + ylim (-4, 4) + xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + mastyle

pcaLS <- grid.arrange(mypca2, pca.load2, nrow = 1)
ggsave("loadScorPCA_EXP2.jpg", plot = pcaLS,
       path = "./results",
       scale = 1, width = 25, height = 20, units = "cm",
       dpi = 600)
################################################################################
# VOLCANO PLOTS ################################################################
res1 <- read.csv(file = "./results/CTR_M3_36S_CLU_KOvsCTR_M3_36S_WT.csv")
res2 <- read.csv(file = "./results/KOLFC2_CLU_KOvsKOLFC2_WT.csv")
res3 <- read.csv(file = "./results/KOLFC2_WTvsCTR_M3_36S_WT.csv")
res4 <- read.csv(file = "./results/KOLFC2_CLU_KOvsCTR_M3_36S_CLU_KO.csv")
################################################################################
top10.res1 <- top_n(res1, -10, pvalue)[, 9]
volcRes1 <- EnhancedVolcano(res1,
                            lab = res1$symbol,
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            selectLab = top10.res1,
                            #xlim = c(-10,10),
                            #ylim = c(0,40),
                            xlab = bquote(~Log[2]~ 'fold-change'),
                            pCutoff = 10e-6,
                            FCcutoff = 1.5,
                            title = "CLU-KO vs WT",
                            subtitle = "geno::CTR.M3.36S",
                            #caption = "FC cutoff, 1.5; p-value cutoff, 10e-6",
                            #pointSize = 4.0,
                            #labSize = 4.0,
                            colAlpha = 1,
                            legendPosition = 'right',
                            legendLabSize = 12,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 0.2,
                            colConnectors = 'grey30',
                            boxedLabels = TRUE,
                            labSize = 4,
                            labCol = 'black',
                            labFace = 'bold') +
    theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none") + labs(caption = NULL)

# ggsave("volcano_CTR_M3_36S_CLU_KOvsCTR_M3_36S_WT.jpg", plot = volcRes1, device = NULL,
#        path = "./results",
#        scale = 1, width = 20, height = 20, units = "cm",
#        dpi = 600)
###
top10.res2 <- top_n(res2, -10, pvalue)[, 9]
volcRes2 <- EnhancedVolcano(res2,
                            lab = res2$symbol,
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            selectLab = top10.res2,
                           # xlim = c(-6,6),
                            xlab = bquote(~Log[2]~ 'fold-change'),
                            pCutoff = 10e-6,
                            FCcutoff = 1.5,
                            title = "CLU-KO vs WT",
                            subtitle = "geno::KOLFC2",
                            #caption = "FC cutoff, 1.333; p-value cutoff, 10e-4",
                            #pointSize = 4.0,
                            #labSize = 4.0,
                            colAlpha = 1,
                            legendPosition = 'right',
                            legendLabSize = 12,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 0.2,
                            colConnectors = 'grey30',
                            boxedLabels = TRUE,
                            labSize = 4,
                            labCol = 'black',
                            labFace = 'bold') +
                theme(plot.title = element_text(hjust = 0.5),
                 plot.subtitle = element_text(hjust = 0.5),
                 legend.position = "none") + labs(caption = NULL)

# ggsave("volcano_KOLFC2_CLU_KOvsKOLFC2_WT.jpg", plot = volcRes2, device = NULL,
#        path = "./results",
#        scale = 1, width = 20, height = 20, units = "cm",
#        dpi = 600)

###
top10.res3 <- top_n(res3, -10, pvalue)[, 9]
volcRes3 <- EnhancedVolcano(res3,
                            lab = res3$symbol,
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            selectLab = top10.res3,
                            #xlim = c(-6,6),
                            xlab = bquote(~Log[2]~ 'fold-change'),
                            pCutoff = 10e-6,
                            FCcutoff = 1.5,
                            title = "KOLFC2 vs CTR.M3.36S",
                            subtitle = "cond::WT",
                            #caption = "FC cutoff, 1.333; p-value cutoff, 10e-4",
                            #pointSize = 4.0,
                            #labSize = 4.0,
                            colAlpha = 1,
                            legendPosition = 'right',
                            legendLabSize = 12,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 0.2,
                            colConnectors = 'grey30',
                            boxedLabels = TRUE,
                            labSize = 4,
                            labCol = 'black',
                            labFace = 'bold') +
    theme(plot.title = element_text(hjust = 0.5),
                 plot.subtitle = element_text(hjust = 0.5),
                 legend.position = "none") + labs(caption = NULL)

# ggsave("volcano_KOLFC2_WTvsCTR_M3_36S_WT.jpg", plot = volcRes3, device = NULL,
#        path = "./results",
#        scale = 0.8, width = 20, height = 20, units = "cm",
#        dpi = 600)
###
top10.res4 <- top_n(res4, -10, pvalue)[, 9]
volcRes4a <- EnhancedVolcano(res4,
                            lab = res4$symbol,
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            selectLab = top10.res4,
                            #xlim = c(-6,6),
                            xlab = bquote(~Log[2]~ 'fold-change'),
                            pCutoff = 10e-6,
                            FCcutoff = 1.5,
                            title = "KOLFC2 vs CTR.M3.36S",
                            subtitle = "cond::CLU-KO",
                            #caption = "FC cutoff, 1.333; p-value cutoff, 10e-4",
                            #pointSize = 4.0,
                            #labSize = 4.0,
                            colAlpha = 1,
                            legendPosition = 'right',
                            legendLabSize = 12,
                            legendIconSize = 4.0,
                            drawConnectors = TRUE,
                            widthConnectors = 0.2,
                            colConnectors = 'grey30',
                            boxedLabels = TRUE,
                            labSize = 4,
                            labCol = 'black',
                            labFace = 'bold') + theme(plot.title = element_text(hjust = 0.5),
                 plot.subtitle = element_text(hjust = 0.5),
                 legend.position = "none") + labs(caption = NULL)

####### extract legend ############################
volcRes4 <- EnhancedVolcano(res4,
                             lab = res4$symbol,
                             x = 'log2FoldChange',
                             y = 'pvalue',
                             selectLab = top10.res4,
                             #xlim = c(-6,6),
                             xlab = bquote(~Log[2]~ 'fold-change'),
                             pCutoff = 10e-6,
                             FCcutoff = 1.5,
                             title = "KOLFC2 vs CTR.M3.36S",
                             subtitle = "cond::CLU-KO",
                             #caption = "FC cutoff, 1.333; p-value cutoff, 10e-4",
                             #pointSize = 4.0,
                             #labSize = 4.0,
                             colAlpha = 1,
                             legendPosition = 'bottom',
                             legendLabSize = 12,
                             legendIconSize = 4.0,
                             drawConnectors = TRUE,
                             widthConnectors = 0.2,
                             colConnectors = 'grey30',
                             boxedLabels = TRUE,
                             labSize = 4,
                             labCol = 'black',
                             labFace = 'bold') + theme(plot.title = element_text(hjust = 0.5),
                                                       plot.subtitle = element_text(hjust = 0.5))
########################################################

# ggsave("volcano_CTR_M3_36S_CRISPRi.CLUvsCRISPRi.Control.jpg", plot = volcRes4, device = NULL,
#        path = "./results",
#        scale = 1, width = 20, height = 20, units = "cm",
#        dpi = 600)
# ###
# ggsave("volcano_m1.jpg", plot = multiplot(volcRes1, volcRes2, cols = 1), device = NULL,
#        path = "F:/data/RNAseq/P190850/analysis/redo",
#        scale = 0.8, width = 20, height = 20, units = "cm",
#        dpi = 600)
#
# ggsave("volcano_m2.jpg", plot = multiplot(volcRes3, volcRes4, cols = 1), device = NULL,
#        path = "F:/data/RNAseq/P190850/analysis/redo",
#        scale = 0.8, width = 20, height = 20, units = "cm",
#        dpi = 600)

library(cowplot)
legend <- get_legend(volcRes4)
volAll <- grid.arrange(volcRes1, volcRes2, volcRes3, volcRes4a, ncol = 2, bottom = "pCutoff = 10e-6 & FCcutoff = 1.5")
vol.All <- plot_grid(volAll, legend, ncol = 1, rel_heights = c(1,0.1))
vol.All
ggsave("volcano_all_exp2.png", plot = vol.All, device = NULL,
path = "./results",
scale = 1, width = 25, height = 30, units = "cm",
dpi = 600)
################################################################################
