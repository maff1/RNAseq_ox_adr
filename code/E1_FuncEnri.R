#!/usr/bin/env Rscript
###############################################################################
## RNAseq workflow - PART E1::Functional Enrichment Analysis ##################
###############################################################################
#BiocManager::install("clusterProfiler")
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

keytypes(org.Hs.eg.db)

res1$entrez <- mapIds(x = org.Hs.eg.db,
                      keys = res1$ensembl,
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")
colnames(res1)
res1.geneList <- res1[, c("entrez","log2FoldChange")]
res1.geneList <- res1.geneList[!duplicated(res1.geneList[, c(1)]),]
res1.geneList <- res1.geneList[!is.na(res1.geneList[, c(1)]),]
res1.geneList <- res1.geneList[!is.na(res1.geneList[, c(2)]),]
row.names(res1.geneList) <- res1.geneList[, c(1)]

res1.geneList <- res1.geneList[order(res1.geneList$log2FoldChange, decreasing = T),]
res1.geneList <- res1.geneList[, c(2), drop=F]

# Extract the foldchanges
foldchanges.res1 <- res1.geneList$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges.res1) <- row.names(res1.geneList)
## Sort fold changes in decreasing order
foldchanges.res1 <- sort(foldchanges.res1, decreasing = TRUE)

options(digits=3)
## GSEA using gene sets from KEGG pathways
res1.gseaKEGG <- gseKEGG(geneList = foldchanges.res1, # ordered named vector of fold changes (Entrez IDs are the associated names)
                         organism = "hsa", # supported organisms listed below
                         nPerm = 1000, # default number permutations
                         minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                         pvalueCutoff = 0.05, # padj cutoff value
                         verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results1 <- res1.gseaKEGG@result

gk1 <- gseKEGG(foldchanges.res1, organism="hsa", nPerm=10000)
head(gk1)
saveRDS(gk1, "./results/gsea_res1.rds")

## dotplot for GSEA result
plot1 <- dotplot(gk1, showCategory=30) +
    facet_grid(.~ifelse(NES < 0, 'suppressed', 'activated'))
plot2 <- gseaplot2(gk1, geneSetID = 1:4)

gkx1 <- setReadable(gk1, 'org.Hs.eg.db', 'ENTREZID')
#cnetplot(gkx, foldChange=foldchanges.res2, showCategory = 1)
plot3 <- emapplot(gk1,layout="kk")
plot4 <- upsetplot(gk1)

gsea.res2.plot <- grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, top = "GSEA geno::KOLFC2 CLU-KO vs WT")
# now add the title
title <- ggdraw() +
    draw_label(
        "GSEA geno::KOLFC2 CLU-KO vs WT",
        fontface = 'bold')

gsea.res2.plot2 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, ncol = 2, labels=LETTERS[1:4])
gsea.res2.plot3 <- cowplot::plot_grid(title, gsea.res2.plot3, ncol = 1, rel_heights = c(0.1, 1))
ggsave("GSEA_KOLFC2_CLUKOvsWT_v2.png", plot = gsea.res2.plot3, device = NULL,
       path = "./results",
       scale = 1, width = 45, height = 30, units = "cm",
       dpi = 600)
###############################################################




## add gene entrez ids
res2$entrez <- mapIds(x = org.Hs.eg.db,
                            keys = res2$ensembl,
                            column = "ENTREZID",
                            keytype = "ENSEMBL",
                            multiVals = "first")
colnames(res2)
res2.geneList <- res2[, c("entrez","log2FoldChange")]
res2.geneList <- res2.geneList[!duplicated(res2.geneList[, c(1)]),]
res2.geneList <- res2.geneList[!is.na(res2.geneList[, c(1)]),]
res2.geneList <- res2.geneList[!is.na(res2.geneList[, c(2)]),]
row.names(res2.geneList) <- res2.geneList[, c(1)]

res2.geneList <- res2.geneList[order(res2.geneList$log2FoldChange, decreasing = T),]
res2.geneList <- res2.geneList[, c(2), drop=F]

# Extract the foldchanges
foldchanges.res2 <- res2.geneList$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges.res2) <- row.names(res2.geneList)
## Sort fold changes in decreasing order
foldchanges.res2 <- sort(foldchanges.res2, decreasing = TRUE)

options(digits=3)
## GSEA using gene sets from KEGG pathways
res2.gseaKEGG <- gseKEGG(geneList = foldchanges.res2, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    nPerm = 1000, # default number permutations
                    minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- res2.gseaKEGG@result

gk <- gseKEGG(foldchanges.res2, organism="hsa", nPerm=10000)
head(gk)
saveRDS(gk, "./results/gsea_res2.rds")

## dotplot for GSEA result
plot1 <- dotplot(gk, showCategory=30) +
    facet_grid(.~ifelse(NES < 0, 'suppressed', 'activated'))
plot2 <- gseaplot2(gk, geneSetID = 1:4)

gkx <- setReadable(gk, 'org.Hs.eg.db', 'ENTREZID')
#cnetplot(gkx, foldChange=foldchanges.res2, showCategory = 1)
plot3 <- emapplot(gk,layout="kk")
plot4 <- upsetplot(gk)

gsea.res2.plot <- grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, top = "GSEA geno::KOLFC2 CLU-KO vs WT")
# now add the title
title <- ggdraw() +
    draw_label(
        "GSEA geno::KOLFC2 CLU-KO vs WT",
        fontface = 'bold')

gsea.res2.plot2 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, ncol = 2, labels=LETTERS[1:4])
gsea.res2.plot3 <- cowplot::plot_grid(title, gsea.res2.plot3, ncol = 1, rel_heights = c(0.1, 1))
ggsave("GSEA_KOLFC2_CLUKOvsWT_v2.png", plot = gsea.res2.plot3, device = NULL,
       path = "./results",
       scale = 1, width = 45, height = 30, units = "cm",
       dpi = 600)
###############################################################
## add gene entrez ids
res3$entrez <- mapIds(x = org.Hs.eg.db,
                      keys = res3$ensembl,
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")
colnames(res3)
res3.geneList <- res3[, c("entrez","log2FoldChange")]
res3.geneList <- res3.geneList[!duplicated(res3.geneList[, c(1)]),]
res3.geneList <- res3.geneList[!is.na(res3.geneList[, c(1)]),]
res3.geneList <- res3.geneList[!is.na(res3.geneList[, c(2)]),]
row.names(res3.geneList) <- res3.geneList[, c(1)]

res3.geneList <- res3.geneList[order(res3.geneList$log2FoldChange, decreasing = T),]
res3.geneList <- res3.geneList[, c(2), drop=F]

# Extract the foldchanges
foldchanges.res3 <- res3.geneList$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges.res3) <- row.names(res3.geneList)
## Sort fold changes in decreasing order
foldchanges.res3 <- sort(foldchanges.res3, decreasing = TRUE)

options(digits=3)
## GSEA using gene sets from KEGG pathways
res3.gseaKEGG <- gseKEGG(geneList = foldchanges.res3, # ordered named vector of fold changes (Entrez IDs are the associated names)
                         organism = "hsa", # supported organisms listed below
                         nPerm = 1000, # default number permutations
                         minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                         pvalueCutoff = 0.05, # padj cutoff value
                         verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results3 <- res3.gseaKEGG@result

gk3 <- gseKEGG(foldchanges.res3, organism="hsa", nPerm=10000)
head(gk3)
saveRDS(gk3, "./results/gsea_res3.rds")

## dotplot for GSEA result
plot1 <- dotplot(gk3, showCategory=30) +
    facet_grid(.~ifelse(NES < 0, 'suppressed', 'activated'))
plot2 <- gseaplot2(gk3, geneSetID = 1:4)

gkx3 <- setReadable(gk3, 'org.Hs.eg.db', 'ENTREZID')
#cnetplot(gkx, foldChange=foldchanges.res2, showCategory = 1)
plot3 <- emapplot(gk3,layout="kk")
plot4 <- upsetplot(gk3)

#gsea.res3.plot <- grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, top = "GSEA cond::WT KOLFC2 vs CTR.M3.36S")
# now add the title
title <- ggdraw() +
    draw_label(
        "GSEA cond::WT KOLFC2 vs CTR.M3.36S",
        fontface = 'bold')

gsea.res3.plot2 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, ncol = 2, labels=LETTERS[1:4])
gsea.res3.plot3 <- cowplot::plot_grid(title, gsea.res3.plot2, ncol = 1, rel_heights = c(0.1, 1))
ggsave("KOLFC2_WTvsCTR_M3_36S_WT_exp2.png", plot = gsea.res3.plot3, device = NULL,
       path = "./results",
       scale = 1, width = 45, height = 30, units = "cm",
       dpi = 600)
###########################################################################
## add gene entrez ids
res4$entrez <- mapIds(x = org.Hs.eg.db,
                      keys = res4$ensembl,
                      column = "ENTREZID",
                      keytype = "ENSEMBL",
                      multiVals = "first")
colnames(res4)
res4.geneList <- res4[, c("entrez","log2FoldChange")]
res4.geneList <- res4.geneList[!duplicated(res4.geneList[, c(1)]),]
res4.geneList <- res4.geneList[!is.na(res4.geneList[, c(1)]),]
res4.geneList <- res4.geneList[!is.na(res4.geneList[, c(2)]),]
row.names(res4.geneList) <- res4.geneList[, c(1)]

res4.geneList <- res4.geneList[order(res4.geneList$log2FoldChange, decreasing = T),]
res4.geneList <- res4.geneList[, c(2), drop=F]

# Extract the foldchanges
foldchanges.res4 <- res4.geneList$log2FoldChange
## Name each fold change with the corresponding Entrez ID
names(foldchanges.res4) <- row.names(res4.geneList)
## Sort fold changes in decreasing order
foldchanges.res4 <- sort(foldchanges.res4, decreasing = TRUE)

options(digits=3)
## GSEA using gene sets from KEGG pathways
res4.gseaKEGG <- gseKEGG(geneList = foldchanges.res4, # ordered named vector of fold changes (Entrez IDs are the associated names)
                         organism = "hsa", # supported organisms listed below
                         nPerm = 1000, # default number permutations
                         minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                         pvalueCutoff = 0.05, # padj cutoff value
                         verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results4 <- res4.gseaKEGG@result

gk4 <- gseKEGG(foldchanges.res4, organism="hsa", nPerm=10000)
head(gk4)
saveRDS(gk4, "./results/gsea_res4.rds")

## dotplot for GSEA result
plot1 <- dotplot(gk4, showCategory=30) +
    facet_grid(.~ifelse(NES < 0, 'suppressed', 'activated'))
plot2 <- gseaplot2(gk4, geneSetID = 1:4)

gkx4 <- setReadable(gk4, 'org.Hs.eg.db', 'ENTREZID')
#cnetplot(gkx, foldChange=foldchanges.res2, showCategory = 1)
plot3 <- emapplot(gk4,layout="kk")
plot4 <- upsetplot(gk4)

#gsea.res3.plot <- grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, top = "GSEA cond::WT KOLFC2 vs CTR.M3.36S")
# now add the title
title <- ggdraw() +
    draw_label(
        "GSEA cond::CLU-KO KOLFC2 vs CTR.M3.36S",
        fontface = 'bold')

gsea.res4.plot2 <- cowplot::plot_grid(plot1, plot2, plot3, plot4, ncol = 2, labels=LETTERS[1:4])
gsea.res4.plot3 <- cowplot::plot_grid(title, gsea.res4.plot2, ncol = 1, rel_heights = c(0.1, 1))
ggsave("KOLFC2_CLU_KOvsCTR_M3_36S_CLU_KO_exp2.png", plot = gsea.res4.plot3, device = NULL,
       path = "./results",
       scale = 1, width = 45, height = 30, units = "cm",
       dpi = 600)
###################################
#BiocManager::install("pathview")
library("pathview")
hsa04151.res3 <- pathview(gene.data  = foldchanges.res3,
                     pathway.id = "hsa04151",
                     species    = "hsa",
                     limit      = list(gene=10, cpd=1))

hsa04621.res3 <- pathview(gene.data  = foldchanges.res3,
                          pathway.id = "hsa04621",
                          species    = "hsa",
                          limit      = list(gene=10, cpd=1))

hsa04510.res3 <- pathview(gene.data  = foldchanges.res3,
                          pathway.id = "hsa04510",
                          species    = "hsa",
                          limit      = list(gene=10, cpd=1))

hsa04670.res3 <- pathview(gene.data  = foldchanges.res3,
                          pathway.id = "hsa04670",
                          species    = "hsa",
                          limit      = list(gene=10, cpd=1))

hsa04210.res3 <- pathview(gene.data  = foldchanges.res3,
                          pathway.id = "hsa04210",
                          species    = "hsa",
                          limit      = list(gene=10, cpd=1))
################################################################################

