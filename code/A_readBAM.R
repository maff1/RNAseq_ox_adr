#!/usr/bin/env Rscript
################################################################################
## RNAseq workflow - PART A:: read BAM and summarise overlaps ##################
################################################################################
packages <- c("readr", "flextable", "R.utils", "GenomicFeatures", "Rsamtools",
              "GenomicAlignments", "BiocManager")
ipak <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE,
                         repos = "http://cran.r-project.org/")
    sapply(pkg, require, character.only = TRUE)
}
ipak(packages)
################################################################################
library(readr)
library(flextable)
library(R.utils)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
################################################################################
setwd("F:/data/PROJECTS/rnaseq/adria")
pathFiles <- "./data/200819_A00711_0245_AHTJLTDMXX/bam"
pathRef <- "./data/200717_A00711_0228_AHTFFFDMXX/referenceGenome"
################################################################################
paths <- file.path(pathFiles, "./sampleTable.txt")
sampleTable <- read_delim(file = paths, delim="\t")
sampleTable
################################################################################
filenames <- file.path(pathFiles, paste0(sampleTable$Run, ".bam"))
file.exists(filenames)
################################################################################
bamfiles <- BamFileList(filenames, yieldSize=2000000)
################################################################################
seqinfo(bamfiles[1])
################################################################################
gtffile <- file.path(pathRef,"Homo_sapiens.GRCh37.87.gtf.gz")
txdb <- makeTxDbFromGFF(gunzip(gtffile, remove=FALSE), format="gtf")
ebg <- exonsBy(txdb, by="gene")
ebg
################################################################################
fls <- list.files(pathFiles, pattern="bam$", full=TRUE)
bamLst <- BamFileList(fls, yieldSize=100000)
################################################################################
sequ <- summarizeOverlaps(ebg,
                         bamLst,
                         mode="Union",
                         singleEnd=FALSE,
                         ignore.strand=TRUE,
                         fragments=TRUE)
saveRDS(sequ, file = "./results/200819_A00711_0245_AHTJLTDMXX_seq.rds")
################################################################################
