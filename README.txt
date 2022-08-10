################################################################################
## README ######################################################################
############################################################## MAFF V.1.0 ######
# TITLE: bulk RNA-seq. analysis P200110
# AUTHOR: Marco Fernandes
# EMAIL: marco(dot)fernandes(at)psych(dot)ox(dot)ac(dot)uk
# DATE: 24 JULY 2020

# PROJECT: P200110
# ORGANISM: Homo sapiens (hsa)
# GENOME: GRCh37.EBVB95-8wt.ERCC
# INSTRUMENT: Illumina ???? ????
# TYPE: RNA-Seq PolyA
# LIBRARY: paired-end | 487/20_MPX
# READ TYPE: two reads and two lanes
# STRAND: no strand specific
# DATA.SOURCE: ftp://gishriwy:darva-vehen-27@bsg-ftp.well.ox.ac.uk/200717_A00711_0228_AHTFFFDMXX
# ORIGINAL.FILES: '.bam', '.bai', '.fastq'

# CONSTRUCTION.PROTOCOL: ??? (to be filled)
# SPIKE-INS: EBVB95-8wt.ERCC https://www.biostars.org/p/339187

# BACKGROUND: sequenced iPSC-neurons' RNA, from cells overexpressing CLU
(CRISPRa) and their respective control cells, and from cells with silencing
of CLU (CRISPRi) and their respective control cells. Each condition done in
triplicate (12 samples) and the whole experiment done in two different
genotypes (total of 24 samples).

# COMPARISONS:
-CLU-CRISPRa vs its control
-CLU-CRISPRi vs its control
-CLU-CRISPRa vs CLU-CRISPRi
-Variability between replicates (CRISPRa-control and CRISPRi-Control are pure
technical replicates, but CRISPRa-CLU and CRISPRi-CLU replicates use different
sgRNAs for modifying CLU's expression)
-Results from KOLFC2 vs CTR_M3_36S genotypes
