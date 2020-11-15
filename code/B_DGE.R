#!/usr/bin/env Rscript
################################################################################
## RNAseq workflow - PART B:: DGE ##############################################
################################################################################
packages <- c("DESeq2", "biomaRt")
ipak <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE,
                         repos = "http://cran.r-project.org/")
    sapply(pkg, require, character.only = TRUE)
}
ipak(packages)
BiocManager::install("BiocParallel")
BiocManager::install("biomaRt")
BiocManager::install(c("org.Hs.eg.db","ensembldb"))
################################################################################
library(DESeq2)
library(biomaRt)
library(ensembldb)
library(stringr)
library(BiocParallel)
library(readr)
################################################################################
#register(MulticoreParam(32))
register(bpstart(SnowParam(4)))   # active snow cluster for the session
################################################################################
if (!exists("sequ")|!exists("sampleTable")) {
    sequ <- readRDS("./results/200819_A00711_0245_AHTJLTDMXX_seq.rds")
} else if (!exists("sampleTable")) {
    sampleTable <- read_delim(list.files(".", recursive = TRUE,
                                         pattern="sampleTable.txt"), delim="\t")
} else {print("summarizeOverlaps and sampleTable already loaded")
}
################################################################################
colnames(sequ) <- gsub("\\.bam", "", colnames(sequ))
seIdx <- match(colnames(sequ), sampleTable$Run)
colData(sequ) <- cbind(Run.1 = colData(sequ), sampleTable[seIdx, ])
################################################################################
sequ$group <- factor(paste0(sequ$Genotype, "_",sequ$Condition))
sequ$Sample <- gsub("-", ".", sequ$Sample)
sequ$Condition <- gsub("-", ".", sequ$Condition)
sequ$group <- gsub("-", ".", sequ$group)
#################################################################################
# assay(sequ)
# colData(sequ)
# rowRanges(sequ)
##################################################################################
dds <- DESeqDataSet(sequ, design = as.formula(~ group))
dds.shr <- DESeq(dds, betaPrior=TRUE, parallel=TRUE)
dds.un.shr <- DESeq(dds, betaPrior=FALSE, parallel=TRUE)

par(mfrow=c(2,1))
plotMA(dds.shr, ylim=c(-10,10), alpha = 0.05, main = "shrunken results")
plotMA(dds.un.shr, ylim=c(-10,10), alpha = 0.05, main = "unshrunken results")
########################
resultsNames(dds.un.shr)
########################
res1 <- results(dds.un.shr, contrast = c("group", "CTR_M3_36S_CLU.KO", "CTR_M3_36S_WT"), parallel=TRUE)
res2 <- results(dds.un.shr, contrast = c("group", "KOLFC2_CLU.KO", "KOLFC2_WT"), parallel=TRUE)
res3 <- results(dds.un.shr, contrast = c("group", "KOLFC2_WT", "CTR_M3_36S_WT"), parallel=TRUE)
res4 <- results(dds.un.shr, contrast = c("group", "KOLFC2_CLU.KO", "CTR_M3_36S_CLU.KO"), parallel=TRUE)
################################################################################
# res1 <- results(dds, contrast = c("group", "KOLFC2_CRISPRa.CLU", "KOLFC2_CRISPRa.Control"))
# res2 <- results(dds, contrast = c("group", "KOLFC2_CRISPRi.CLU", "KOLFC2_CRISPRi.Control"))
# res3 <- results(dds, contrast = c("group", "CTR_M3_36S_CRISPRa.CLU", "CTR_M3_36S_CRISPRa.Control"))
# res4 <- results(dds, contrast = c("group", "CTR_M3_36S_CRISPRi.CLU", "CTR_M3_36S_CRISPRi.Control"))
# #
# res5 <- results(dds, contrast = c("group", "KOLFC2_CRISPRa.CLU", "CTR_M3_36S_CRISPRa.CLU"))
# res6 <- results(dds, contrast = c("group", "KOLFC2_CRISPRi.CLU", "CTR_M3_36S_CRISPRi.CLU"))
# res7 <- results(dds, contrast = c("group", "KOLFC2_CRISPRa.Control", "CTR_M3_36S_CRISPRa.Control"))
# res8 <- results(dds, contrast = c("group", "KOLFC2_CRISPRi.Control", "CTR_M3_36S_CRISPRi.Control"))
################################################################################
saveRDS(dds.un.shr, file = "./results/dExp_200819_A00711_0245_AHTJLTDMXX.rds")
dds.un.shr <- readRDS(file = "./results/dExp_200819_A00711_0245_AHTJLTDMXX.rds")
################################################################################
res1$ensembl <- sapply( strsplit( rownames(res1), split="\\+" ), "[", 1 )
res2$ensembl <- sapply( strsplit( rownames(res2), split="\\+" ), "[", 1 )
res3$ensembl <- sapply( strsplit( rownames(res3), split="\\+" ), "[", 1 )
res4$ensembl <- sapply( strsplit( rownames(res4), split="\\+" ), "[", 1 )
################################################################################
# res1$ensembl <- sapply( strsplit( rownames(res1), split="\\+" ), "[", 1 )
# res2$ensembl <- sapply( strsplit( rownames(res2), split="\\+" ), "[", 1 )
# res3$ensembl <- sapply( strsplit( rownames(res3), split="\\+" ), "[", 1 )
# res4$ensembl <- sapply( strsplit( rownames(res4), split="\\+" ), "[", 1 )
# res5$ensembl <- sapply( strsplit( rownames(res5), split="\\+" ), "[", 1 )
# res6$ensembl <- sapply( strsplit( rownames(res6), split="\\+" ), "[", 1 )
# res7$ensembl <- sapply( strsplit( rownames(res6), split="\\+" ), "[", 1 )
# res8$ensembl <- sapply( strsplit( rownames(res6), split="\\+" ), "[", 1 )
################################################################################
# ensembl <- useEnsembl(biomart = "ensembl",
#                       dataset = "hsapiens_gene_ensembl",
#                       mirror = "useast")
# ensembl = useMart("ensembl",
#                   host="useast.ensembl.org",
#                   ensemblRedirect = FALSE, dataset = "hsapiens_gene_ensembl")
#
# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# genemap <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
#                   filters = "ensembl_gene_id",
#                   values = res1$ensembl,
#                   mart = ensembl)

#########################
gtffile <- file.path(pathRef,"Homo_sapiens.GRCh37.87.gtf.gz")
## generate the SQLite database file
DB <- ensDbFromGtf(gtf = gtffile)
## load the DB file directly
EDB <- EnsDb(DB)
keytypes(EDB)
###############################
res1$symbol <- mapIds(x = EDB,
                      keys = res1$ensembl,
                      column = "SYMBOL",
                      keytype = "GENEID",
                      multiVals = "first")

res1$genebiotype <- mapIds(x = EDB,
                      keys = res1$ensembl,
                      column = "GENEBIOTYPE",
                      keytype = "GENEID",
                      multiVals = "first")
###############################
res2$symbol <- mapIds(x = EDB,
                      keys = res2$ensembl,
                      column = "SYMBOL",
                      keytype = "GENEID",
                      multiVals = "first")

res2$genebiotype <- mapIds(x = EDB,
                           keys = res2$ensembl,
                           column = "GENEBIOTYPE",
                           keytype = "GENEID",
                           multiVals = "first")
###############################
res3$symbol <- mapIds(x = EDB,
                      keys = res3$ensembl,
                      column = "SYMBOL",
                      keytype = "GENEID",
                      multiVals = "first")

res3$genebiotype <- mapIds(x = EDB,
                           keys = res3$ensembl,
                           column = "GENEBIOTYPE",
                           keytype = "GENEID",
                           multiVals = "first")
###############################
res4$symbol <- mapIds(x = EDB,
                      keys = res4$ensembl,
                      column = "SYMBOL",
                      keytype = "GENEID",
                      multiVals = "first")

res4$genebiotype <- mapIds(x = EDB,
                           keys = res4$ensembl,
                           column = "GENEBIOTYPE",
                           keytype = "GENEID",
                           multiVals = "first")
#####################################
# t.res1 <- res1[1:10, ]
# t.res2 <- res2[1:10, ]
#
# mapS <- function(dat) {
# dat$symbol <- mapIds(x = EDB,
#            keys = dat$ensembl,
#            column = "SYMBOL",
#            keytype = "GENEID",
#            multiVals = "first")
# }
# mls <- list(t.res1, t.res2)
# lapply(mls, mapS)
################################################################################
# idx1 <- match(res1$ensembl, genemap$ensembl_gene_id)
# res1$entrez <- genemap$entrezgene[idx1]
# res1$hgnc_symbol <- genemap$hgnc_symbol[idx1]
# ###
# idx2 <- match(res2$ensembl, genemap$ensembl_gene_id)
# res2$entrez <- genemap$entrezgene[idx2]
# res2$hgnc_symbol <- genemap$hgnc_symbol[idx2]
# ###
# idx3 <- match(res3$ensembl, genemap$ensembl_gene_id)
# res3$entrez <- genemap$entrezgene[idx3]
# res3$hgnc_symbol <- genemap$hgnc_symbol[idx3]
# ###
# idx4 <- match(res4$ensembl, genemap$ensembl_gene_id)
# res4$entrez <- genemap$entrezgene[idx4]
# res4$hgnc_symbol <- genemap$hgnc_symbol[idx4]
###
# idx5 <- match(res5$ensembl, genemap$ensembl_gene_id)
# res5$entrez <- genemap$entrezgene[idx5]
# res5$hgnc_symbol <- genemap$hgnc_symbol[idx5]
# ###
# idx6 <- match(res6$ensembl, genemap$ensembl_gene_id)
# res6$entrez <- genemap$entrezgene[idx6]
# res6$hgnc_symbol <- genemap$hgnc_symbol[idx6]
# ###
# idx7 <- match(res7$ensembl, genemap$ensembl_gene_id)
# res7$entrez <- genemap$entrezgene[idx7]
# res7$hgnc_symbol <- genemap$hgnc_symbol[idx7]
# ###
# idx8 <- match(res8$ensembl, genemap$ensembl_gene_id)
# res8$entrez <- genemap$entrezgene[idx8]
# res8$hgnc_symbol <- genemap$hgnc_symbol[idx8]
################################################################################
# t.res1 <- res1[1:10, ]
# t.res2 <- res2[1:10, ]
# t.res3 <- res2[1:10, ]
# t.res4 <- res2[1:10, ]
#
# l.res <- list(t.res1,t.res2,t.res3,t.res4)
# names(l.res) <- c("1", "2", "3", "4")
# names(l.res)
# sapply(names(l.res),
#        function (x) write.csv(l.res[[x]][order(l.res[[x]]$padj,
#                                                decreasing = FALSE, na.last = TRUE), ],
#                               file=paste0("./results/", x, ".csv")))
################################################################################
l.res <- list(res1,res2,res3,res4)
names(l.res) <- c("CTR_M3_36S_CLU_KOvsCTR_M3_36S_WT",
                  "KOLFC2_CLU_KOvsKOLFC2_WT",
                  "KOLFC2_WTvsCTR_M3_36S_WT",
                  "KOLFC2_CLU_KOvsCTR_M3_36S_CLU_KO")
names(l.res)
sapply(names(l.res),
       function (x) write.csv(l.res[[x]][order(l.res[[x]]$padj,
                                               decreasing = FALSE, na.last = TRUE), ],
                              file=paste0("./results/", x, ".csv")))

################################################################################
# res1$Length <- str_count(res1$hgnc_symbol)
# res1$hgnc_symbol <- ifelse(is.na(res1$hgnc_symbol) | res1$Length == 0,
#                            res1$ensembl, res1$hgnc_symbol)
#
# write.csv(as.data.frame(res1[order(res1$padj, decreasing = FALSE, na.last = TRUE), ]),
#           file = "./results/KOLFC2_CRISPRa.CLUvsCRISPRa.Control.csv")
# ###
# res2$Length <- str_count(res2$hgnc_symbol)
# res2$hgnc_symbol <- ifelse(is.na(res2$hgnc_symbol) | res2$Length == 0,
#                            res2$ensembl, res2$hgnc_symbol)
#
# write.csv(as.data.frame(res2[order(res2$padj, decreasing = FALSE, na.last = TRUE), ]),
#           file = "./results/KOLFC2_CRISPRi.CLUvsCRISPRi.Control.csv")
#
# ###
# res3$Length <- str_count(res3$hgnc_symbol)
# res3$hgnc_symbol <- ifelse(is.na(res3$hgnc_symbol) | res3$Length == 0,
#                            res3$ensembl, res3$hgnc_symbol)
#
# write.csv(as.data.frame(res3[order(res3$padj, decreasing = FALSE, na.last = TRUE), ]),
#           file = "./results/CTR_M3_36S_CRISPRa.CLUvsCRISPRa.Control.csv")
#
# ###
# res4$Length <- str_count(res4$hgnc_symbol)
# res4$hgnc_symbol <- ifelse(is.na(res4$hgnc_symbol) | res4$Length == 0,
#                            res4$ensembl, res4$hgnc_symbol)
#
# write.csv(as.data.frame(res4[order(res4$padj, decreasing = FALSE, na.last = TRUE), ]),
#           file = "./results/CTR_M3_36S_CRISPRi.CLUvsCRISPRi.Control.csv")
# ###
# res5$Length <- str_count(res5$hgnc_symbol)
# res5$hgnc_symbol <- ifelse(is.na(res5$hgnc_symbol) | res5$Length == 0,
#                            res5$ensembl, res5$hgnc_symbol)
#
# write.csv(as.data.frame(res5[order(res5$padj, decreasing = FALSE, na.last = TRUE), ]),
#           file = "./results/KOLFC2_CRISPRa.CLUvsCTR_M3_36S_CRISPRa.CLU.csv")
# ###
# res6$Length <- str_count(res6$hgnc_symbol)
# res6$hgnc_symbol <- ifelse(is.na(res6$hgnc_symbol) | res6$Length == 0,
#                            res6$ensembl, res6$hgnc_symbol)
#
# write.csv(as.data.frame(res6[order(res6$padj, decreasing = FALSE, na.last = TRUE), ]),
#           file = "./results/KOLFC2_CRISPRi.CLUvsCTR_M3_36S_CRISPRi.CLU.csv")
# ###
# res7$Length <- str_count(res7$hgnc_symbol)
# res7$hgnc_symbol <- ifelse(is.na(res7$hgnc_symbol) | res7$Length == 0,
#                            res7$ensembl, res7$hgnc_symbol)
#
# write.csv(as.data.frame(res7[order(res7$padj, decreasing = FALSE, na.last = TRUE), ]),
#           file = "./results/KOLFC2_CRISPRa.ControlvsCTR_M3_36S_CRISPRa.Control.csv")
# ###
# res8$Length <- str_count(res8$hgnc_symbol)
# res8$hgnc_symbol <- ifelse(is.na(res8$hgnc_symbol) | res8$Length == 0,
#                            res8$ensembl, res8$hgnc_symbol)
#
# write.csv(as.data.frame(res8[order(res8$padj, decreasing = FALSE, na.last = TRUE), ]),
#           file = "./results/KOLFC2_CRISPRi.ControlvsCTR_M3_36S_CRISPRi.Control.csv")
#################################################################################

## TODO #########################################################################
## Extract condition effects by adding interactions (linear combinations) ######
# https://support.bioconductor.org/p/65087/
# https://support.bioconductor.org/p/85823/
# dds$group <- factor(paste0(dds$genotype, dds$condition))
# design(dds) <- ~ group
# dds <- DESeq(dds)
# resultsNames(dds)
# results(dds, contrast=c("group", "IB", "IA"))
################################################################################
