#!/usr/bin/env Rscript
################################################################################
## RNAseq workflow - PART C:: diagnostic.QC.post.DEG ###########################
################################################################################
packages <- c("DESeq2", "dplyr", "reshape2", "ggplot2")
ipak <- function(pkg) {
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE,
                         repos = "http://cran.r-project.org/")
    sapply(pkg, require, character.only = TRUE)
}
ipak(packages)
BiocManager::install("vsn")
################################################################################
library(DESeq2)
library(vsn)
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
################################################################################
fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
}
## VARIANCE EFFECTS ############################################################
vsd <- vst(dds.un.shr, blind=FALSE)
rld <- rlog(dds.un.shr, blind=FALSE)

# this gives log2(n + 1)
ntd <- normTransform(dds.un.shr, f = log2, pc = 1)
# this gives log2(n + 100)
ntd100 <- normTransform(dds.un.shr, f = log2, pc = 100)

p.ntd <- meanSdPlot(assay(ntd))
p.ntd100 <- meanSdPlot(assay(ntd100))
p.vsd <- meanSdPlot(assay(vsd))
p.rld <- meanSdPlot(assay(rld))

## VARIANCE EFFECTS PLOTS ########################################################
p1 <- p.ntd$gg + ggtitle("ntd, log(counts+1)") + scale_fill_gradient(low = "blue", high = "darkred") + scale_y_continuous(limits = c(0, 4)) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title=element_text(size=12,face="bold"),
          legend.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5, linetype = "solid"))

p2 <- p.vsd$gg + ggtitle("vsd") + scale_fill_gradient(low = "blue", high = "darkred") + scale_y_continuous(limits = c(0, 4)) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title=element_text(size=12,face="bold"),
          legend.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5, linetype = "solid"))

p3 <- p.rld$gg + ggtitle("rld") + scale_fill_gradient(low = "blue", high = "darkred") + scale_y_continuous(limits = c(0, 4)) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title=element_text(size=12,face="bold"),
          legend.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5, linetype = "solid"))

p4 <- p.ntd100$gg + ggtitle("ntd, log(counts+100)") + scale_fill_gradient(low = "blue", high = "darkred") + scale_y_continuous(limits = c(0, 4)) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title=element_text(size=12,face="bold"),
          legend.title = element_text(face = "bold"),
          panel.background = element_rect(fill = "lightblue",
                                          colour = "lightblue",
                                          size = 0.5, linetype = "solid"))

(layout_matrix <- matrix(c(1, 1, 2, 2, 4, 3, 3, 4), nrow = 2, byrow = TRUE))
eda1 <- grid.arrange(p1, p3, p2, layout_matrix = layout_matrix,
             top = textGrob("EDA: variance-stabilizing transformations",
                            gp=gpar(fontsize=18,font=2)))

eda2 <- grid.arrange(p1, p4, p3, p2, ncol = 2, nrow = 2,
                     top = textGrob("EDA: variance-stabilizing transformations",
                                    gp=gpar(fontsize=18,font=2)))

ggsave("./results/eda1_exp2.png", plot = eda1, width = 7, height = 7, units = "in", dpi = 300)
ggsave("./results/eda2_exp2.png", plot = eda2, width = 7, height = 7, units = "in", dpi = 300)
################################################################################
## DIAGNOSTIC PLOTS ############################################################
plotDispEsts2 <- function(dds){

    as.data.frame(mcols(dds)) %>%
        dplyr::select(baseMean, dispGeneEst, dispFit, dispersion) %>%
        melt(id.vars="baseMean") %>%
        dplyr::filter(baseMean>0) %>%
        ggplot(aes(x=baseMean, y=value, colour=variable)) +
        geom_point(size=0.1) +
        scale_x_log10() +
        scale_y_log10(labels=fancy_scientific) +
        theme_bw() +
        ylab("Dispersion") +
        xlab("BaseMean") +
        scale_colour_manual(
            values=c("Black", "#e41a1c", "#377eb8"),
            breaks=c("dispGeneEst", "dispFit", "dispersion"),
            labels=c("Estimate", "Fit", "Final"),
            name=""
        ) +
        guides(colour = guide_legend(override.aes = list(size=2)))
}
#dds.un.shr
plotDispEsts2(dds.un.shr)
################################################################################
png(filename = "dispEst.png", units="in", width=5, height=5, res=300)
plotDispEsts(dds.un.shr)
dev.off()
#plotDispEsts2(dds)

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds.un.shr)[["cooks"]]), range=0, las=2)


resGA <- results(dds.un.shr, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds.un.shr, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds.un.shr, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds.un.shr, lfcThreshold=.5, altHypothesis="less")

pdf("plotMA_EXP2.pdf", width=10,height=10)
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-10,10)
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim, xlab="mean expression", ylab="log2fc"); drawLines()
plotMA(resL, ylim=ylim); drawLines()
mtext("Wald tests for log2FC thresholdings", outer = TRUE, cex = 1.5)
dev.off()
###########

# resApeT <- lfcShrink(dds.un.shr, coef=2, type="apeglm", lfcThreshold=1)
# plotMA(resApeT, ylim=c(-3,3), cex=.8)
# abline(h=c(-1,1), col="dodgerblue", lwd=2)
################################################################################
