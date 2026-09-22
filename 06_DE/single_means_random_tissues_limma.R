#!/usr/bin/env Rscript

# R script to perform DE analysis with DREAM

# turn off scientific notation
options(scipen = 999)

# import libraries
#library(variancePartition)
library(edgeR)
#library(BiocParallel)
library(tibble)
library(rcartocolor)
library(ggVennDiagram)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random_tissues"
dir.create(workingDir)
setwd(workingDir)

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/MackLab/metabolic_adaptation/counted/counts_merged_samples.csv", row.names="gene")

# import grouping factor
factors <- read.csv(file="/Users/bamflappy/MackLab/input_data/experimental_design.csv", row.names="sample")

# set input FDR cutoff
cutFDR <- 0.05

# set LFC cut
#cutFC <- 1.2
cutLFC <- log2(1.2)

# trim the data table to remove lines with counting statistics (htseq)
removeList <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
countsTable <- inputTable[!row.names(inputTable) %in% removeList,]

# convert the grouping data into factors 
targets <- as.data.frame(lapply(factors, as.factor))
rownames(targets) <- rownames(factors)

# specify the reference level
#targets$treatment <- relevel(targets$treatment, ref = "Ad_lib")
#targets$mouseline <- relevel(targets$mouseline, ref = "SARA")
#targets$tissue <- relevel(targets$tissue, ref = "BAT")
#targets$sex <- relevel(targets$sex, ref = "M")

# sort samples alphabetically
countsTable <- countsTable[, order(colnames(countsTable))]
targets <- targets[order(rownames(targets)), ]

# check sample naming
setdiff(colnames(countsTable), rownames(targets))
setdiff(rownames(targets), colnames(countsTable))

# setup a design matrix
sample_group <- factor(paste(targets$mouseline, targets$tissue, targets$sex, targets$treatment, sep="."))

# specify the reference level
sample_group <- relevel(sample_group, ref = "SARA.BAT.Male.Ad_lib")

# there are at least 2 samples per individual
smallestGroupSize <- 2

# keep only rows that have a count of at least 10 for a minimal number of samples
keep <- rowSums(cpm(countsTable) >= 0.1) >= smallestGroupSize
list <- DGEList(countsTable[keep, ])

# percentage of genes lost from filtering
(nrow(countsTable) - nrow(list)) / nrow(countsTable) * 100

# use TMM normalization to eliminate composition biases between libraries
#dge <- calcNormFactors(list)
dge <- normLibSizes(list)
# write normalized counts to file
normList <- cpm(dge, normalized.lib.sizes=TRUE)
# add gene row name tag
normList <- as_tibble(normList, rownames = "gene")
write.table(normList, file="normalizedCounts.csv", sep=",", row.names=FALSE, quote=FALSE)

# write log transformed normalized counts to file
normListLog <- cpm(dge, normalized.lib.sizes=TRUE, log=TRUE)
normListLog <- as_tibble(normListLog, rownames = "gene")
write.table(normListLog, file="normalizedCounts_logTransformed.csv", sep=",", row.names=FALSE, quote=FALSE)


## limma analysis

# experimental design as a one-way layout, where one coefficient is assigned to each group
design <- model.matrix(~ 0 + sample_group)
colnames(design) <- levels(sample_group)

# view the design
design

# apply duplicateCorrelation is two rounds
vobj_tmp <- voom(dge, design, plot = FALSE)
dupcor <- duplicateCorrelation(vobj_tmp, design, block=targets$individual)

# view correlation between repeated measurements
dupcor$consensus.correlation

# run voom again to compute more accurate precision weights
vobj <- voom(dge, design, plot = FALSE, block = targets$individual, correlation = dupcor$consensus)

# estimate linear mixed model with a single variance component and fit the model for each gene
dupcor <- duplicateCorrelation(vobj, design, block = targets$individual)

# this step uses only the genome-wide average for the random effect
fitDupCor <- lmFit(vobj, design, block = targets$individual, correlation = dupcor$consensus)

# view the factors
colnames(design)

# create lists of factors
factor_names <- colnames(design)
factor_Ad <- factor_names[grepl("Ad_lib", factor_names)]
factor_Food <- factor_names[grepl("Food_Restriction", factor_names)]
factor_F <- factor_names[grepl("Female", factor_names)]
factor_M <- factor_names[grepl("Male", factor_names)]
factor_M_Ad <- factor_M[grepl("Ad_lib", factor_M)]
factor_M_Food <- factor_M[grepl("Food_Restriction", factor_M)]
factor_SARA <- factor_names[grepl("SARA", factor_names)]
factor_SARA_Ad <- factor_SARA[grepl("Ad_lib", factor_SARA)]
factor_SARA_Food <- factor_SARA[grepl("Food_Restriction", factor_SARA)]
factor_SARA_F <- factor_SARA[grepl("Female", factor_SARA)]
factor_SARA_M <- factor_SARA[grepl("Male", factor_SARA)]
factor_SARA_M_Ad <- factor_SARA_M[grepl("Ad_lib", factor_SARA_M)]
factor_SARA_M_Food <- factor_SARA_M[grepl("Food_Restriction", factor_SARA_M)]
factor_SARA_BAT <- factor_SARA[grepl("BAT", factor_SARA)]
factor_SARA_BAT_Ad <- factor_SARA_BAT[grepl("Ad_lib", factor_SARA_BAT)]
factor_SARA_BAT_M_Ad <- factor_SARA_BAT_Ad[grepl("Male", factor_SARA_BAT_Ad)]
factor_SARA_BAT_Food <- factor_SARA_BAT[grepl("Food_Restriction", factor_SARA_BAT)]
factor_SARA_WAT <- factor_SARA[grepl("WAT", factor_SARA)]
factor_SARA_WAT_Ad <- factor_SARA_WAT[grepl("Ad_lib", factor_SARA_WAT)]
factor_SARA_WAT_M_Ad <- factor_SARA_WAT_Ad[grepl("Male", factor_SARA_WAT_Ad)]
factor_SARA_WAT_Food <- factor_SARA_WAT[grepl("Food_Restriction", factor_SARA_WAT)]
factor_SARA_HYP <- factor_SARA[grepl("HYP", factor_SARA)]
factor_SARA_HYP_Ad <- factor_SARA_HYP[grepl("Ad_lib", factor_SARA_HYP)]
factor_SARA_HYP_M_Ad <- factor_SARA_HYP_Ad[grepl("Male", factor_SARA_HYP_Ad)]
factor_SARA_HYP_Food <- factor_SARA_HYP[grepl("Food_Restriction", factor_SARA_HYP)]
factor_SARA_LIV <- factor_SARA[grepl("LIV", factor_SARA)]
factor_SARA_LIV_Ad <- factor_SARA_LIV[grepl("Ad_lib", factor_SARA_LIV)]
factor_SARA_LIV_M_Ad <- factor_SARA_LIV_Ad[grepl("Male", factor_SARA_LIV_Ad)]
factor_SARA_LIV_Food <- factor_SARA_LIV[grepl("Food_Restriction", factor_SARA_LIV)]
factor_MANF <- factor_names[grepl("MANF", factor_names)]
factor_MANF_Ad <- factor_MANF[grepl("Ad_lib", factor_MANF)]
factor_MANF_Food <- factor_MANF[grepl("Food_Restriction", factor_MANF)]
factor_MANF_F <- factor_MANF[grepl("Female", factor_MANF)]
factor_MANF_M <- factor_MANF[grepl("Male", factor_MANF)]
factor_MANF_M_Ad <- factor_MANF_M[grepl("Ad_lib", factor_MANF_M)]
factor_MANF_M_Food <- factor_MANF_M[grepl("Food_Restriction", factor_MANF_M)]
factor_MANF_BAT <- factor_MANF[grepl("BAT", factor_MANF)]
factor_MANF_BAT_Ad <- factor_MANF_BAT[grepl("Ad_lib", factor_MANF_BAT)]
factor_MANF_BAT_M_Ad <- factor_MANF_BAT_Ad[grepl("Male", factor_MANF_BAT_Ad)]
factor_MANF_BAT_Food <- factor_MANF_BAT[grepl("Food_Restriction", factor_MANF_BAT)]
factor_MANF_WAT <- factor_MANF[grepl("WAT", factor_MANF)]
factor_MANF_WAT_Ad <- factor_MANF_WAT[grepl("Ad_lib", factor_MANF_WAT)]
factor_MANF_WAT_M_Ad <- factor_MANF_WAT_Ad[grepl("Male", factor_MANF_WAT_Ad)]
factor_MANF_WAT_Food <- factor_MANF_WAT[grepl("Food_Restriction", factor_MANF_WAT)]
factor_MANF_HYP <- factor_MANF[grepl("HYP", factor_MANF)]
factor_MANF_HYP_Ad <- factor_MANF_HYP[grepl("Ad_lib", factor_MANF_HYP)]
factor_MANF_HYP_M_Ad <- factor_MANF_HYP_Ad[grepl("Male", factor_MANF_HYP_Ad)]
factor_MANF_HYP_Food <- factor_MANF_HYP[grepl("Food_Restriction", factor_MANF_HYP)]
factor_MANF_LIV <- factor_MANF[grepl("LIV", factor_MANF)]
factor_MANF_LIV_Ad <- factor_MANF_LIV[grepl("Ad_lib", factor_MANF_LIV)]
factor_MANF_LIV_M_Ad <- factor_MANF_LIV_Ad[grepl("Male", factor_MANF_LIV_Ad)]
factor_MANF_LIV_Food <- factor_MANF_LIV[grepl("Food_Restriction", factor_MANF_LIV)]

# setup contrasts
contrasts <- makeContrasts(
  SARA_BAT_M_AdvsFood=
    eval(parse(text=paste(factor_SARA_BAT_M_Ad, collapse = "+")))/length(factor_SARA_BAT_M_Ad) -
    eval(parse(text=paste(factor_SARA_BAT_Food, collapse = "+")))/length(factor_SARA_BAT_Food),
  SARA_WAT_M_AdvsFood=
    eval(parse(text=paste(factor_SARA_WAT_M_Ad, collapse = "+")))/length(factor_SARA_WAT_M_Ad) -
    eval(parse(text=paste(factor_SARA_WAT_Food, collapse = "+")))/length(factor_SARA_WAT_Food),
  SARA_HYP_M_AdvsFood=
    eval(parse(text=paste(factor_SARA_HYP_M_Ad, collapse = "+")))/length(factor_SARA_HYP_M_Ad) -
    eval(parse(text=paste(factor_SARA_HYP_Food, collapse = "+")))/length(factor_SARA_HYP_Food),
  SARA_LIV_M_AdvsFood=
    eval(parse(text=paste(factor_SARA_LIV_M_Ad, collapse = "+")))/length(factor_SARA_LIV_M_Ad) -
    eval(parse(text=paste(factor_SARA_LIV_Food, collapse = "+")))/length(factor_SARA_LIV_Food),
  MANF_BAT_AdvsFood=
    eval(parse(text=paste(factor_MANF_BAT_Ad, collapse = "+")))/length(factor_MANF_BAT_Ad) -
    eval(parse(text=paste(factor_MANF_BAT_Food, collapse = "+")))/length(factor_MANF_BAT_Food),
  MANF_WAT_AdvsFood=
    eval(parse(text=paste(factor_MANF_WAT_Ad, collapse = "+")))/length(factor_MANF_WAT_Ad) -
    eval(parse(text=paste(factor_MANF_WAT_Food, collapse = "+")))/length(factor_MANF_WAT_Food),
  MANF_HYP_AdvsFood=
    eval(parse(text=paste(factor_MANF_HYP_Ad, collapse = "+")))/length(factor_MANF_HYP_Ad) -
    eval(parse(text=paste(factor_MANF_HYP_Food, collapse = "+")))/length(factor_MANF_HYP_Food),
  MANF_LIV_AdvsFood=
    eval(parse(text=paste(factor_MANF_LIV_Ad, collapse = "+")))/length(factor_MANF_LIV_Ad) -
    eval(parse(text=paste(factor_MANF_LIV_Food, collapse = "+")))/length(factor_MANF_LIV_Food),
  MANF_BAT_M_AdvsFood=
    eval(parse(text=paste(factor_MANF_BAT_M_Ad, collapse = "+")))/length(factor_MANF_BAT_M_Ad) -
    eval(parse(text=paste(factor_MANF_BAT_Food, collapse = "+")))/length(factor_MANF_BAT_Food),
  MANF_WAT_M_AdvsFood=
    eval(parse(text=paste(factor_MANF_WAT_M_Ad, collapse = "+")))/length(factor_MANF_WAT_M_Ad) -
    eval(parse(text=paste(factor_MANF_WAT_Food, collapse = "+")))/length(factor_MANF_WAT_Food),
  MANF_HYP_M_AdvsFood=
    eval(parse(text=paste(factor_MANF_HYP_M_Ad, collapse = "+")))/length(factor_MANF_HYP_M_Ad) -
    eval(parse(text=paste(factor_MANF_HYP_Food, collapse = "+")))/length(factor_MANF_HYP_Food),
  MANF_LIV_M_AdvsFood=
    eval(parse(text=paste(factor_MANF_LIV_M_Ad, collapse = "+")))/length(factor_MANF_LIV_M_Ad) -
    eval(parse(text=paste(factor_MANF_LIV_Food, collapse = "+")))/length(factor_MANF_LIV_Food),
  SARAvsMANF_BAT_M_AdvsFood=
    (eval(parse(text=paste(factor_SARA_BAT_M_Ad, collapse = "+")))/length(factor_SARA_BAT_M_Ad) -
    eval(parse(text=paste(factor_SARA_BAT_Food, collapse = "+")))/length(factor_SARA_BAT_Food)) -
    (eval(parse(text=paste(factor_MANF_BAT_M_Ad, collapse = "+")))/length(factor_MANF_BAT_M_Ad) -
       eval(parse(text=paste(factor_MANF_BAT_Food, collapse = "+")))/length(factor_MANF_BAT_Food)),
  SARAvsMANF_WAT_M_AdvsFood=
    (eval(parse(text=paste(factor_SARA_WAT_M_Ad, collapse = "+")))/length(factor_SARA_WAT_M_Ad) -
    eval(parse(text=paste(factor_SARA_WAT_Food, collapse = "+")))/length(factor_SARA_WAT_Food)) -
    (eval(parse(text=paste(factor_MANF_WAT_M_Ad, collapse = "+")))/length(factor_MANF_WAT_M_Ad) -
       eval(parse(text=paste(factor_MANF_WAT_Food, collapse = "+")))/length(factor_MANF_WAT_Food)),
  SARAvsMANF_HYP_M_AdvsFood=
    (eval(parse(text=paste(factor_SARA_HYP_M_Ad, collapse = "+")))/length(factor_SARA_HYP_M_Ad) -
    eval(parse(text=paste(factor_SARA_HYP_Food, collapse = "+")))/length(factor_SARA_HYP_Food)) -
    (eval(parse(text=paste(factor_MANF_HYP_M_Ad, collapse = "+")))/length(factor_MANF_HYP_M_Ad) -
       eval(parse(text=paste(factor_MANF_HYP_Food, collapse = "+")))/length(factor_MANF_HYP_Food)),
  SARAvsMANF_LIV_M_AdvsFood=
    (eval(parse(text=paste(factor_SARA_LIV_M_Ad, collapse = "+")))/length(factor_SARA_LIV_M_Ad) -
    eval(parse(text=paste(factor_SARA_LIV_Food, collapse = "+")))/length(factor_SARA_LIV_Food)) -
    (eval(parse(text=paste(factor_MANF_LIV_M_Ad, collapse = "+")))/length(factor_MANF_LIV_M_Ad) -
       eval(parse(text=paste(factor_MANF_LIV_Food, collapse = "+")))/length(factor_MANF_LIV_Food)),
  levels=colnames(design)) 

# view contrasts
#plotContrasts(contrasts)
#contrasts

# fit the contrasts
fitCont <- contrasts.fit(fitDupCor, contrasts)

# fit Empirical Bayes for moderated t-statistics
fitCont <- eBayes(fitCont)

# get names of available coefficients and contrasts for testing
colnames(fitCont)

# treatment: Food_Restriction vs Ad_lib (12-hr food restriction vs. unlimited food, male mice only)
# tissue: LIV (liver), HYP (hypothalamus), BAT (brown adipose), WAT (white adipose)
# mouseline: MANF vs SARA (skinny vs fat)
# sex: M vs F (unlimited food)

# How do male SARA tissues respond to 12-hr food restriction?
## SARA_BAT_M_AdvsFood
## SARA_WAT_M_AdvsFood
## SARA_HYP_M_AdvsFood
## SARA_LIV_M_AdvsFood
# How do male MANF tissues respond to 12-hr food restriction?
## MANF_BAT_M_AdvsFood
## MANF_WAT_M_AdvsFood
## MANF_HYP_M_AdvsFood
## MANF_LIV_M_AdvsFood
# How do SARA M respond differently than MANF M to 12-hr food restriction? 
## SARAvsMANF_BAT_M_AdvsFood
## SARAvsMANF_WAT_M_AdvsFood
## SARAvsMANF_HYP_M_AdvsFood
## SARAvsMANF_LIV_M_AdvsFood

# get genes < FDR and LFC cutoffs
# treatment among male tissues
treatment_SARA_BAT_M_dge_sig <- topTable(fitCont, coef = "SARA_BAT_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_WAT_M_dge_sig <- topTable(fitCont, coef = "SARA_WAT_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_HYP_M_dge_sig <- topTable(fitCont, coef = "SARA_HYP_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_LIV_M_dge_sig <- topTable(fitCont, coef = "SARA_LIV_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_BAT_M_dge_sig <- topTable(fitCont, coef = "MANF_BAT_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_WAT_M_dge_sig <- topTable(fitCont, coef = "MANF_WAT_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_HYP_M_dge_sig <- topTable(fitCont, coef = "MANF_HYP_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_LIV_M_dge_sig <- topTable(fitCont, coef = "MANF_LIV_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_BAT_M_dge_sig <- topTable(fitCont, coef = "SARAvsMANF_BAT_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_WAT_M_dge_sig <- topTable(fitCont, coef = "SARAvsMANF_WAT_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_HYP_M_dge_sig <- topTable(fitCont, coef = "SARAvsMANF_HYP_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_LIV_M_dge_sig <- topTable(fitCont, coef = "SARAvsMANF_LIV_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)

# check the number of sig DE genes
# treatment among male tissues
nrow(treatment_SARA_BAT_M_dge_sig)
nrow(treatment_SARA_WAT_M_dge_sig)
nrow(treatment_SARA_HYP_M_dge_sig)
nrow(treatment_SARA_LIV_M_dge_sig)
nrow(treatment_MANF_BAT_M_dge_sig)
nrow(treatment_MANF_WAT_M_dge_sig)
nrow(treatment_MANF_HYP_M_dge_sig)
nrow(treatment_MANF_LIV_M_dge_sig)
nrow(treatment_SARA_MANF_BAT_M_dge_sig)
nrow(treatment_SARA_MANF_WAT_M_dge_sig)
nrow(treatment_SARA_MANF_HYP_M_dge_sig)
nrow(treatment_SARA_MANF_LIV_M_dge_sig)

# export tables of sig DE genes
# male tissue SARA: treatment_SARA_BAT_M_dge_sig, treatment_SARA_WAT_M_dge_sig, treatment_SARA_HYP_M_dge_sig, treatment_SARA_LIV_M_dge_sig
treatment_SARA_BAT_M_dge_sig_tbl <- as_tibble(treatment_SARA_BAT_M_dge_sig, rownames = "gene")
treatment_SARA_BAT_M_out_file <- paste("treatment_SARA_BAT_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_BAT_M_out_file <- paste(treatment_SARA_BAT_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_BAT_M_dge_sig_tbl, file=treatment_SARA_BAT_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_WAT_M_dge_sig_tbl <- as_tibble(treatment_SARA_WAT_M_dge_sig, rownames = "gene")
treatment_SARA_WAT_M_out_file <- paste("treatment_SARA_WAT_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_WAT_M_out_file <- paste(treatment_SARA_WAT_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_WAT_M_dge_sig_tbl, file=treatment_SARA_WAT_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_HYP_M_dge_sig_tbl <- as_tibble(treatment_SARA_HYP_M_dge_sig, rownames = "gene")
treatment_SARA_HYP_M_out_file <- paste("treatment_SARA_HYP_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_HYP_M_out_file <- paste(treatment_SARA_HYP_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_HYP_M_dge_sig_tbl, file=treatment_SARA_HYP_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_LIV_M_dge_sig_tbl <- as_tibble(treatment_SARA_LIV_M_dge_sig, rownames = "gene")
treatment_SARA_LIV_M_out_file <- paste("treatment_SARA_LIV_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_LIV_M_out_file <- paste(treatment_SARA_LIV_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_LIV_M_dge_sig_tbl, file=treatment_SARA_LIV_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
# male tissue MANF: treatment_MANF_BAT_M_dge_sig, treatment_MANF_WAT_M_dge_sig, treatment_MANF_HYP_M_dge_sig, treatment_MANF_LIV_M_dge_sig
treatment_MANF_BAT_M_dge_sig_tbl <- as_tibble(treatment_MANF_BAT_M_dge_sig, rownames = "gene")
treatment_MANF_BAT_M_out_file <- paste("treatment_MANF_BAT_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_MANF_BAT_M_out_file <- paste(treatment_MANF_BAT_M_out_file, "csv", sep = ".")
write.table(treatment_MANF_BAT_M_dge_sig_tbl, file=treatment_MANF_BAT_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_WAT_M_dge_sig_tbl <- as_tibble(treatment_MANF_WAT_M_dge_sig, rownames = "gene")
treatment_MANF_WAT_M_out_file <- paste("treatment_MANF_WAT_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_MANF_WAT_M_out_file <- paste(treatment_MANF_WAT_M_out_file, "csv", sep = ".")
write.table(treatment_MANF_WAT_M_dge_sig_tbl, file=treatment_MANF_WAT_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_HYP_M_dge_sig_tbl <- as_tibble(treatment_MANF_HYP_M_dge_sig, rownames = "gene")
treatment_MANF_HYP_M_out_file <- paste("treatment_MANF_HYP_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_MANF_HYP_M_out_file <- paste(treatment_MANF_HYP_M_out_file, "csv", sep = ".")
write.table(treatment_MANF_HYP_M_dge_sig_tbl, file=treatment_MANF_HYP_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_LIV_M_dge_sig_tbl <- as_tibble(treatment_MANF_LIV_M_dge_sig, rownames = "gene")
treatment_MANF_LIV_M_out_file <- paste("treatment_MANF_LIV_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_MANF_LIV_M_out_file <- paste(treatment_MANF_LIV_M_out_file, "csv", sep = ".")
write.table(treatment_MANF_LIV_M_dge_sig_tbl, file=treatment_MANF_LIV_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
# male tissue SARA vs MANF: treatment_SARA_MANF_BAT_M_dge_sig, treatment_SARA_MANF_WAT_M_dge_sig, treatment_SARA_MANF_HYP_M_dge_sig, treatment_SARA_MANF_LIV_M_dge_sig
treatment_SARA_MANF_BAT_M_dge_sig_tbl <- as_tibble(treatment_SARA_MANF_BAT_M_dge_sig, rownames = "gene")
treatment_SARA_MANF_BAT_M_out_file <- paste("treatment_SARA_MANF_BAT_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_MANF_BAT_M_out_file <- paste(treatment_SARA_MANF_BAT_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_MANF_BAT_M_dge_sig_tbl, file=treatment_SARA_MANF_BAT_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_WAT_M_dge_sig_tbl <- as_tibble(treatment_SARA_MANF_WAT_M_dge_sig, rownames = "gene")
treatment_SARA_MANF_WAT_M_out_file <- paste("treatment_SARA_MANF_WAT_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_MANF_WAT_M_out_file <- paste(treatment_SARA_MANF_WAT_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_MANF_WAT_M_dge_sig_tbl, file=treatment_SARA_MANF_WAT_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_HYP_M_dge_sig_tbl <- as_tibble(treatment_SARA_MANF_HYP_M_dge_sig, rownames = "gene")
treatment_SARA_MANF_HYP_M_out_file <- paste("treatment_SARA_MANF_HYP_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_MANF_HYP_M_out_file <- paste(treatment_SARA_MANF_HYP_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_MANF_HYP_M_dge_sig_tbl, file=treatment_SARA_MANF_HYP_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_LIV_M_dge_sig_tbl <- as_tibble(treatment_SARA_MANF_LIV_M_dge_sig, rownames = "gene")
treatment_SARA_MANF_LIV_M_out_file <- paste("treatment_SARA_MANF_LIV_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_MANF_LIV_M_out_file <- paste(treatment_SARA_MANF_LIV_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_MANF_LIV_M_dge_sig_tbl, file=treatment_SARA_MANF_LIV_M_out_file, sep=",", row.names=FALSE, quote=FALSE)

# get all results of hypothesis tests
# male tissue SARA: treatment_SARA_BAT_M_dge_sig, treatment_SARA_WAT_M_dge_sig, treatment_SARA_HYP_M_dge_sig, treatment_SARA_LIV_M_dge_sig
treatment_SARA_BAT_M_dge_results <- topTable(fitCont, coef = "SARA_BAT_M_AdvsFood", number = nrow(dge))
treatment_SARA_WAT_M_dge_results <- topTable(fitCont, coef = "SARA_WAT_M_AdvsFood", number = nrow(dge))
treatment_SARA_HYP_M_dge_results <- topTable(fitCont, coef = "SARA_HYP_M_AdvsFood", number = nrow(dge))
treatment_SARA_LIV_M_dge_results <- topTable(fitCont, coef = "SARA_LIV_M_AdvsFood", number = nrow(dge))
# male tissue MANF: treatment_MANF_BAT_M_dge_sig, treatment_MANF_WAT_M_dge_sig, treatment_MANF_HYP_M_dge_sig, treatment_MANF_LIV_M_dge_sig
treatment_MANF_BAT_M_dge_results <- topTable(fitCont, coef = "MANF_BAT_M_AdvsFood", number = nrow(dge))
treatment_MANF_WAT_M_dge_results <- topTable(fitCont, coef = "MANF_WAT_M_AdvsFood", number = nrow(dge))
treatment_MANF_HYP_M_dge_results <- topTable(fitCont, coef = "MANF_HYP_M_AdvsFood", number = nrow(dge))
treatment_MANF_LIV_M_dge_results <- topTable(fitCont, coef = "MANF_LIV_M_AdvsFood", number = nrow(dge))
# male tissue SARA vs MANF: treatment_SARA_MANF_BAT_M_dge_sig, treatment_SARA_MANF_WAT_M_dge_sig, treatment_SARA_MANF_HYP_M_dge_sig, treatment_SARA_MANF_LIV_M_dge_sig
treatment_SARA_MANF_BAT_M_dge_results <- topTable(fitCont, coef = "SARAvsMANF_BAT_M_AdvsFood", number = nrow(dge))
treatment_SARA_MANF_WAT_M_dge_results <- topTable(fitCont, coef = "SARAvsMANF_WAT_M_AdvsFood", number = nrow(dge))
treatment_SARA_MANF_HYP_M_dge_results <- topTable(fitCont, coef = "SARAvsMANF_HYP_M_AdvsFood", number = nrow(dge))
treatment_SARA_MANF_LIV_M_dge_results <- topTable(fitCont, coef = "SARAvsMANF_LIV_M_AdvsFood", number = nrow(dge))

# export table of DE genes
# male tissue SARA: treatment_SARA_BAT_M_dge_sig, treatment_SARA_WAT_M_dge_sig, , treatment_SARA_LIV_M_dge_sig
treatment_SARA_BAT_M_dge_results_tbl <- as_tibble(treatment_SARA_BAT_M_dge_results, rownames = "gene")
write.table(treatment_SARA_BAT_M_dge_results_tbl, file="treatment_SARA_BAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_WAT_M_dge_results_tbl <- as_tibble(treatment_SARA_WAT_M_dge_results, rownames = "gene")
write.table(treatment_SARA_WAT_M_dge_results_tbl, file="treatment_SARA_WAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_HYP_M_dge_results_tbl <- as_tibble(treatment_SARA_HYP_M_dge_results, rownames = "gene")
write.table(treatment_SARA_HYP_M_dge_results_tbl, file="treatment_SARA_HYP_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_LIV_M_dge_results_tbl <- as_tibble(treatment_SARA_LIV_M_dge_results, rownames = "gene")
write.table(treatment_SARA_LIV_M_dge_results_tbl, file="treatment_SARA_LIV_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
# male tissue MANF: treatment_MANF_BAT_M_dge_sig, treatment_MANF_WAT_M_dge_sig, , treatment_MANF_LIV_M_dge_sig
treatment_MANF_BAT_M_dge_results_tbl <- as_tibble(treatment_MANF_BAT_M_dge_results, rownames = "gene")
write.table(treatment_MANF_BAT_M_dge_results_tbl, file="treatment_MANF_BAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_WAT_M_dge_results_tbl <- as_tibble(treatment_MANF_WAT_M_dge_results, rownames = "gene")
write.table(treatment_MANF_WAT_M_dge_results_tbl, file="treatment_MANF_WAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_HYP_M_dge_results_tbl <- as_tibble(treatment_MANF_HYP_M_dge_results, rownames = "gene")
write.table(treatment_MANF_HYP_M_dge_results_tbl, file="treatment_MANF_HYP_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_LIV_M_dge_results_tbl <- as_tibble(treatment_MANF_LIV_M_dge_results, rownames = "gene")
write.table(treatment_MANF_LIV_M_dge_results_tbl, file="treatment_MANF_LIV_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
# male tissue SARA vs MANF: treatment_SARA_MANF_BAT_M_dge_sig, treatment_SARA_MANF_WAT_M_dge_sig, , treatment_SARA_MANF_LIV_M_dge_sig
treatment_SARA_MANF_BAT_M_dge_results_tbl <- as_tibble(treatment_SARA_MANF_BAT_M_dge_results, rownames = "gene")
write.table(treatment_SARA_MANF_BAT_M_dge_results_tbl, file="treatment_SARA_MANF_BAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_WAT_M_dge_results_tbl <- as_tibble(treatment_SARA_MANF_WAT_M_dge_results, rownames = "gene")
write.table(treatment_SARA_MANF_WAT_M_dge_results_tbl, file="treatment_SARA_MANF_WAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_HYP_M_dge_results_tbl <- as_tibble(treatment_SARA_MANF_HYP_M_dge_results, rownames = "gene")
write.table(treatment_SARA_MANF_HYP_M_dge_results_tbl, file="treatment_SARA_MANF_HYP_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_LIV_M_dge_results_tbl <- as_tibble(treatment_SARA_MANF_LIV_M_dge_results, rownames = "gene")
write.table(treatment_SARA_MANF_LIV_M_dge_results_tbl, file="treatment_SARA_MANF_LIV_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)

# subset counts table by DE gene set
# SARA
treatment_SARA_BAT_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_SARA_BAT_M_dge_sig_tbl$gene
treatment_SARA_BAT_M_dge_sig_tbl_logcounts <- normListLog[treatment_SARA_BAT_M_dge_sig_tbl.keep, ]
treatment_SARA_WAT_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_SARA_WAT_M_dge_sig_tbl$gene
treatment_SARA_WAT_M_dge_sig_tbl_logcounts <- normListLog[treatment_SARA_WAT_M_dge_sig_tbl.keep, ]
treatment_SARA_HYP_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_SARA_HYP_M_dge_sig_tbl$gene
treatment_SARA_HYP_M_dge_sig_tbl_logcounts <- normListLog[treatment_SARA_HYP_M_dge_sig_tbl.keep, ]
treatment_SARA_LIV_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_SARA_LIV_M_dge_sig_tbl$gene
treatment_SARA_LIV_M_dge_sig_tbl_logcounts <- normListLog[treatment_SARA_LIV_M_dge_sig_tbl.keep, ]
# MANF
treatment_MANF_BAT_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_MANF_BAT_M_dge_sig_tbl$gene
treatment_MANF_BAT_M_dge_sig_tbl_logcounts <- normListLog[treatment_MANF_BAT_M_dge_sig_tbl.keep, ]
treatment_MANF_WAT_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_MANF_WAT_M_dge_sig_tbl$gene
treatment_MANF_WAT_M_dge_sig_tbl_logcounts <- normListLog[treatment_MANF_WAT_M_dge_sig_tbl.keep, ]
treatment_MANF_HYP_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_MANF_HYP_M_dge_sig_tbl$gene
treatment_MANF_HYP_M_dge_sig_tbl_logcounts <- normListLog[treatment_MANF_HYP_M_dge_sig_tbl.keep, ]
treatment_MANF_LIV_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_MANF_LIV_M_dge_sig_tbl$gene
treatment_MANF_LIV_M_dge_sig_tbl_logcounts <- normListLog[treatment_MANF_LIV_M_dge_sig_tbl.keep, ]
# SARA vs MANF
treatment_SARA_MANF_BAT_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_SARA_MANF_BAT_M_dge_sig_tbl$gene
treatment_SARA_MANF_BAT_M_dge_sig_tbl_logcounts <- normListLog[treatment_SARA_MANF_BAT_M_dge_sig_tbl.keep, ]
treatment_SARA_MANF_WAT_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_SARA_MANF_WAT_M_dge_sig_tbl$gene
treatment_SARA_MANF_WAT_M_dge_sig_tbl_logcounts <- normListLog[treatment_SARA_MANF_WAT_M_dge_sig_tbl.keep, ]
treatment_SARA_MANF_HYP_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_SARA_MANF_HYP_M_dge_sig_tbl$gene
treatment_SARA_MANF_HYP_M_dge_sig_tbl_logcounts <- normListLog[treatment_SARA_MANF_HYP_M_dge_sig_tbl.keep, ]
treatment_SARA_MANF_LIV_M_dge_sig_tbl.keep <- normListLog$gene %in% treatment_SARA_MANF_LIV_M_dge_sig_tbl$gene
treatment_SARA_MANF_LIV_M_dge_sig_tbl_logcounts <- normListLog[treatment_SARA_MANF_LIV_M_dge_sig_tbl.keep, ]

# format for plotting
# SARA
treatment_SARA_BAT_M_names <- treatment_SARA_BAT_M_dge_sig_tbl_logcounts$gene
treatment_SARA_BAT_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_SARA_BAT_M_dge_sig_tbl_logcounts) <- treatment_SARA_BAT_M_names
treatment_SARA_WAT_M_names <- treatment_SARA_WAT_M_dge_sig_tbl_logcounts$gene
treatment_SARA_WAT_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_SARA_WAT_M_dge_sig_tbl_logcounts) <- treatment_SARA_WAT_M_names
treatment_SARA_HYP_M_names <- treatment_SARA_HYP_M_dge_sig_tbl_logcounts$gene
treatment_SARA_HYP_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_SARA_HYP_M_dge_sig_tbl_logcounts) <- treatment_SARA_HYP_M_names
treatment_SARA_LIV_M_names <- treatment_SARA_LIV_M_dge_sig_tbl_logcounts$gene
treatment_SARA_LIV_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_SARA_LIV_M_dge_sig_tbl_logcounts) <- treatment_SARA_LIV_M_names
# MANF
treatment_MANF_BAT_M_names <- treatment_MANF_BAT_M_dge_sig_tbl_logcounts$gene
treatment_MANF_BAT_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_MANF_BAT_M_dge_sig_tbl_logcounts) <- treatment_MANF_BAT_M_names
treatment_MANF_WAT_M_names <- treatment_MANF_WAT_M_dge_sig_tbl_logcounts$gene
treatment_MANF_WAT_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_MANF_WAT_M_dge_sig_tbl_logcounts) <- treatment_MANF_WAT_M_names
treatment_MANF_HYP_M_names <- treatment_MANF_HYP_M_dge_sig_tbl_logcounts$gene
treatment_MANF_HYP_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_MANF_HYP_M_dge_sig_tbl_logcounts) <- treatment_MANF_HYP_M_names
treatment_MANF_LIV_M_names <- treatment_MANF_LIV_M_dge_sig_tbl_logcounts$gene
treatment_MANF_LIV_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_MANF_LIV_M_dge_sig_tbl_logcounts) <- treatment_MANF_LIV_M_names
# SARA vs MANF
treatment_SARA_MANF_BAT_M_names <- treatment_SARA_MANF_BAT_M_dge_sig_tbl_logcounts$gene
treatment_SARA_MANF_BAT_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_SARA_MANF_BAT_M_dge_sig_tbl_logcounts) <- treatment_SARA_MANF_BAT_M_names
treatment_SARA_MANF_WAT_M_names <- treatment_SARA_MANF_WAT_M_dge_sig_tbl_logcounts$gene
treatment_SARA_MANF_WAT_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_SARA_MANF_WAT_M_dge_sig_tbl_logcounts) <- treatment_SARA_MANF_WAT_M_names
treatment_SARA_MANF_HYP_M_names <- treatment_SARA_MANF_HYP_M_dge_sig_tbl_logcounts$gene
treatment_SARA_MANF_HYP_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_SARA_MANF_HYP_M_dge_sig_tbl_logcounts) <- treatment_SARA_MANF_HYP_M_names
treatment_SARA_MANF_LIV_M_names <- treatment_SARA_MANF_LIV_M_dge_sig_tbl_logcounts$gene
treatment_SARA_MANF_LIV_M_dge_sig_tbl_logcounts$gene <- NULL
rownames(treatment_SARA_MANF_LIV_M_dge_sig_tbl_logcounts) <- treatment_SARA_MANF_LIV_M_names

# heatmap of results
# SARA
jpeg("treatment_SARA_BAT_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_SARA_BAT_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("treatment_SARA_WAT_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_SARA_WAT_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("treatment_SARA_HYP_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_SARA_HYP_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("treatment_SARA_LIV_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_SARA_LIV_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
# MANF
jpeg("treatment_MANF_BAT_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_MANF_BAT_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("treatment_MANF_WAT_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_MANF_WAT_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("treatment_MANF_HYP_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_MANF_HYP_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("treatment_MANF_LIV_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_MANF_LIV_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
# SARA vs MANF
jpeg("treatment_SARA_MANF_BAT_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_SARA_MANF_BAT_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("treatment_SARA_MANF_WAT_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_SARA_MANF_WAT_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("treatment_SARA_MANF_HYP_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_SARA_MANF_HYP_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("treatment_SARA_MANF_LIV_M_dge_sig_tbl_logcounts.jpg")
heatmap(as.matrix(treatment_SARA_MANF_LIV_M_dge_sig_tbl_logcounts), margins = c(8, 1), labRow = FALSE)
dev.off()

# sort results by gene name
# SARA
treatment_SARA_BAT_M_dge_results_tbl <- treatment_SARA_BAT_M_dge_results_tbl[order(treatment_SARA_BAT_M_dge_results_tbl$gene), ]
treatment_SARA_WAT_M_dge_results_tbl <- treatment_SARA_WAT_M_dge_results_tbl[order(treatment_SARA_WAT_M_dge_results_tbl$gene), ]
treatment_SARA_HYP_M_dge_results_tbl <- treatment_SARA_HYP_M_dge_results_tbl[order(treatment_SARA_HYP_M_dge_results_tbl$gene), ]
treatment_SARA_LIV_M_dge_results_tbl <- treatment_SARA_LIV_M_dge_results_tbl[order(treatment_SARA_LIV_M_dge_results_tbl$gene), ]
# MANF
treatment_MANF_BAT_M_dge_results_tbl <- treatment_MANF_BAT_M_dge_results_tbl[order(treatment_MANF_BAT_M_dge_results_tbl$gene), ]
treatment_MANF_WAT_M_dge_results_tbl <- treatment_MANF_WAT_M_dge_results_tbl[order(treatment_MANF_WAT_M_dge_results_tbl$gene), ]
treatment_MANF_HYP_M_dge_results_tbl <- treatment_MANF_HYP_M_dge_results_tbl[order(treatment_MANF_HYP_M_dge_results_tbl$gene), ]
treatment_MANF_LIV_M_dge_results_tbl <- treatment_MANF_LIV_M_dge_results_tbl[order(treatment_MANF_LIV_M_dge_results_tbl$gene), ]

# store LFCs
gene_IDs <- treatment_SARA_BAT_M_dge_results_tbl$gene
init_DE <- rep("None", nrow(treatment_SARA_BAT_M_dge_results_tbl))
BAT_logFC <- data.frame(
  gene = gene_IDs,
  SARA_LFC = treatment_SARA_BAT_M_dge_results_tbl$logFC,
  MANF_LFC = treatment_MANF_BAT_M_dge_results_tbl$logFC,
  SARA_DE = init_DE,
  MANF_DE = init_DE,
  Response = init_DE
)
WAT_logFC <- data.frame(
  gene = gene_IDs,
  SARA_LFC = treatment_SARA_WAT_M_dge_results_tbl$logFC,
  MANF_LFC = treatment_MANF_WAT_M_dge_results_tbl$logFC,
  SARA_DE = init_DE,
  MANF_DE = init_DE,
  Response = init_DE
)
HYP_logFC <- data.frame(
  gene = gene_IDs,
  SARA_LFC = treatment_SARA_HYP_M_dge_results_tbl$logFC,
  MANF_LFC = treatment_MANF_HYP_M_dge_results_tbl$logFC,
  SARA_DE = init_DE,
  MANF_DE = init_DE,
  Response = init_DE
)
LIV_logFC <- data.frame(
  gene = gene_IDs,
  SARA_LFC = treatment_SARA_LIV_M_dge_results_tbl$logFC,
  MANF_LFC = treatment_MANF_LIV_M_dge_results_tbl$logFC,
  SARA_DE = init_DE,
  MANF_DE = init_DE,
  Response = init_DE
)
# indicate DE gene set
# SARA
treatment_SARA_BAT_M_dge_sig_tbl.keep <- BAT_logFC$gene %in% treatment_SARA_BAT_M_dge_sig_tbl$gene
BAT_logFC$SARA_DE[treatment_SARA_BAT_M_dge_sig_tbl.keep] <- "SARA"
treatment_SARA_WAT_M_dge_sig_tbl.keep <- WAT_logFC$gene %in% treatment_SARA_WAT_M_dge_sig_tbl$gene
WAT_logFC$SARA_DE[treatment_SARA_WAT_M_dge_sig_tbl.keep] <- "SARA"
treatment_SARA_HYP_M_dge_sig_tbl.keep <- HYP_logFC$gene %in% treatment_SARA_HYP_M_dge_sig_tbl$gene
HYP_logFC$SARA_DE[treatment_SARA_HYP_M_dge_sig_tbl.keep] <- "SARA"
treatment_SARA_LIV_M_dge_sig_tbl.keep <- LIV_logFC$gene %in% treatment_SARA_LIV_M_dge_sig_tbl$gene
LIV_logFC$SARA_DE[treatment_SARA_LIV_M_dge_sig_tbl.keep] <- "SARA"
# MANF
treatment_MANF_BAT_M_dge_sig_tbl.keep <- BAT_logFC$gene %in% treatment_MANF_BAT_M_dge_sig_tbl$gene
BAT_logFC$MANF_DE[treatment_MANF_BAT_M_dge_sig_tbl.keep] <- "MANF"
treatment_MANF_WAT_M_dge_sig_tbl.keep <- WAT_logFC$gene %in% treatment_MANF_WAT_M_dge_sig_tbl$gene
WAT_logFC$MANF_DE[treatment_MANF_WAT_M_dge_sig_tbl.keep] <- "MANF"
treatment_MANF_HYP_M_dge_sig_tbl.keep <- HYP_logFC$gene %in% treatment_MANF_HYP_M_dge_sig_tbl$gene
HYP_logFC$MANF_DE[treatment_MANF_HYP_M_dge_sig_tbl.keep] <- "MANF"
treatment_MANF_LIV_M_dge_sig_tbl.keep <- LIV_logFC$gene %in% treatment_MANF_LIV_M_dge_sig_tbl$gene
LIV_logFC$MANF_DE[treatment_MANF_LIV_M_dge_sig_tbl.keep] <- "MANF"
# both
BAT_logFC$Response <- paste(BAT_logFC$SARA_DE, BAT_logFC$MANF_DE, sep = "_")
WAT_logFC$Response <- paste(WAT_logFC$SARA_DE, WAT_logFC$MANF_DE, sep = "_")
HYP_logFC$Response <- paste(HYP_logFC$SARA_DE, HYP_logFC$MANF_DE, sep = "_")
LIV_logFC$Response <- paste(LIV_logFC$SARA_DE, LIV_logFC$MANF_DE, sep = "_")
BAT_logFC$Response <- gsub("SARA_MANF", "Both", BAT_logFC$Response)
WAT_logFC$Response <- gsub("SARA_MANF", "Both", WAT_logFC$Response)
HYP_logFC$Response <- gsub("SARA_MANF", "Both", HYP_logFC$Response)
LIV_logFC$Response <- gsub("SARA_MANF", "Both", LIV_logFC$Response)
# clean up
BAT_logFC$Response <- gsub("None_None", "None", BAT_logFC$Response)
WAT_logFC$Response <- gsub("None_None", "None", WAT_logFC$Response)
HYP_logFC$Response <- gsub("None_None", "None", HYP_logFC$Response)
LIV_logFC$Response <- gsub("None_None", "None", LIV_logFC$Response)
BAT_logFC$Response <- gsub("SARA_None", "SARA", BAT_logFC$Response)
WAT_logFC$Response <- gsub("SARA_None", "SARA", WAT_logFC$Response)
HYP_logFC$Response <- gsub("SARA_None", "SARA", HYP_logFC$Response)
LIV_logFC$Response <- gsub("SARA_None", "SARA", LIV_logFC$Response)
BAT_logFC$Response <- gsub("None_MANF", "MANF", BAT_logFC$Response)
WAT_logFC$Response <- gsub("None_MANF", "MANF", WAT_logFC$Response)
HYP_logFC$Response <- gsub("None_MANF", "MANF", HYP_logFC$Response)
LIV_logFC$Response <- gsub("None_MANF", "MANF", LIV_logFC$Response)

# order data frames
custom_order <- c("None", "SARA", "MANF", "Both")
BAT_logFC_sorted <- BAT_logFC %>%
  arrange(match(Response, custom_order))
WAT_logFC_sorted <- WAT_logFC %>%
  arrange(match(Response, custom_order))
HYP_logFC_sorted <- HYP_logFC %>%
  arrange(match(Response, custom_order))
LIV_logFC_sorted <- LIV_logFC %>%
  arrange(match(Response, custom_order))

# scatter plots of log2 FC
jpeg("treatment_SARA_MANF_BAT_M_dge_LFC.jpg")
ggplot(BAT_logFC_sorted, aes(x = SARA_LFC, y = MANF_LFC, color = Response)) + 
  geom_point(aes(alpha = Response)) +
  scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "MANF" = 0.75)) + 
  scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[4], "MANF" = plotColors[6])) +
  labs(x = "SARA log2 FC", y = "MANF log2 FC") +
  theme_minimal()
dev.off()
jpeg("treatment_SARA_MANF_WAT_M_dge_LFC.jpg")
ggplot(WAT_logFC_sorted, aes(x = SARA_LFC, y = MANF_LFC, color = Response)) + 
  geom_point(aes(alpha = Response)) +
  scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "MANF" = 0.75, "Both" = 1.0)) + 
  scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[4], "MANF" = plotColors[6], "Both" = plotColors[5])) +
  labs(x = "SARA log2 FC", y = "MANF log2 FC") +
  theme_minimal()
dev.off()
jpeg("treatment_SARA_MANF_HYP_M_dge_LFC.jpg")
ggplot(HYP_logFC_sorted, aes(x = SARA_LFC, y = MANF_LFC, color = Response)) + 
  geom_point(aes(alpha = Response)) +
  scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "MANF" = 0.75, "Both" = 1.0)) + 
  scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[4], "MANF" = plotColors[6], "Both" = plotColors[5])) +
  labs(x = "SARA log2 FC", y = "MANF log2 FC") +
  theme_minimal()
dev.off()
jpeg("treatment_SARA_MANF_LIV_M_dge_LFC.jpg")
ggplot(LIV_logFC_sorted, aes(x = SARA_LFC, y = MANF_LFC, color = Response)) + 
  geom_point(aes(alpha = Response)) +
  scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "MANF" = 0.75, "Both" = 1.0)) + 
  scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[4], "MANF" = plotColors[6], "Both" = plotColors[5])) +
  labs(x = "SARA log2 FC", y = "MANF log2 FC") +
  theme_minimal()
dev.off()

# turn on scientific notation
options(scipen = 0)

# scatter plots of log2 FC with correlations
jpeg("treatment_SARA_MANF_BAT_M_dge_LFC_cor.jpg")
ggplot(BAT_logFC_sorted, aes(x = SARA_LFC, y = MANF_LFC, color = Response)) + 
  geom_point(aes(alpha = Response)) +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  stat_cor(method = "pearson") +
  scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "MANF" = 0.75)) + 
  scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[4], "MANF" = plotColors[6])) +
  labs(x = "SARA log2 FC", y = "MANF log2 FC") +
  theme_minimal()
dev.off()
jpeg("treatment_SARA_MANF_WAT_M_dge_LFC_cor.jpg")
ggplot(WAT_logFC_sorted, aes(x = SARA_LFC, y = MANF_LFC, color = Response)) + 
  geom_point(aes(alpha = Response)) +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  stat_cor(method = "pearson") +
  scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "MANF" = 0.75, "Both" = 1.0)) + 
  scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[4], "MANF" = plotColors[6], "Both" = plotColors[5])) +
  labs(x = "SARA log2 FC", y = "MANF log2 FC") +
  theme_minimal()
dev.off()
jpeg("treatment_SARA_MANF_HYP_M_dge_LFC_cor.jpg")
ggplot(HYP_logFC_sorted, aes(x = SARA_LFC, y = MANF_LFC, color = Response)) + 
  geom_point(aes(alpha = Response)) +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  stat_cor(method = "pearson") +
  scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "MANF" = 0.75, "Both" = 1.0)) + 
  scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[4], "MANF" = plotColors[6], "Both" = plotColors[5])) +
  labs(x = "SARA log2 FC", y = "MANF log2 FC") +
  theme_minimal()
dev.off()
jpeg("treatment_SARA_MANF_LIV_M_dge_LFC_cor.jpg")
ggplot(LIV_logFC_sorted, aes(x = SARA_LFC, y = MANF_LFC, color = Response)) + 
  geom_point(aes(alpha = Response)) +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  stat_cor(method = "pearson") +
  scale_alpha_manual(values = c("None" = 0.2, "SARA" = 0.75, "MANF" = 0.75, "Both" = 1.0)) + 
  scale_color_manual(values = c("None" = "#BEBEBE", "SARA" = plotColors[4], "MANF" = plotColors[6], "Both" = plotColors[5])) +
  labs(x = "SARA log2 FC", y = "MANF log2 FC") +
  theme_minimal()
dev.off()

