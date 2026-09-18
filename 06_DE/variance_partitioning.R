#!/usr/bin/env Rscript

# R script to perform variance partitioning

# turn off scientific notation
options(scipen = 999)

# import libraries
library(edgeR)
library(tibble)
library(dplyr)
library(variancePartition)
library(BiocParallel)

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/variance_partitioning"
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

# specify parallel processing parameters
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~ tissue + treatment + sex + mouseline + mouseline:(1|individual) + tissue:mouseline + tissue:sex + treatment:tissue + treatment:mouseline + mouseline:sex

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, targets, BPPARAM = param)

# check the variance partitioning to identify important variables that should be included as fixed or random effects
formVP <- ~ tissue + treatment + sex + mouseline + tissue:mouseline + tissue:sex + treatment:tissue + treatment:mouseline + mouseline:sex
vp <- fitExtractVarPartModel(vobjDream, formVP, targets)

# violin plot of contribution of each variable to total variance
jpeg("varPart.jpg")
plotVarPart(sortCols(vp))
dev.off()

# The variable to be tested must be a fixed effect
form <- ~ tissue + mouseline + tissue:mouseline + mouseline:(1|individual)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, targets, BPPARAM = param)

# check the variance partitioning to identify important variables that should be included as fixed or random effects
formVP <- ~ tissue + mouseline + tissue:mouseline
vp <- fitExtractVarPartModel(vobjDream, formVP, targets)

# violin plot of contribution of each variable to total variance
jpeg("varPart_tissue_mouseline.jpg")
plotVarPart(sortCols(vp))
dev.off()

# check the variance partitioning to identify important variables that should be included as fixed or random effects
formVP <- ~ (1|tissue) + (1|mouseline) + (1|tissue:mouseline) + (1|individual) + (1|mouseline:individual)
vp <- fitExtractVarPartModel(vobjDream, formVP, targets)

# violin plot of contribution of each variable to total variance
jpeg("varPart_random_tissue_mouseline.jpg")
plotVarPart(sortCols(vp))
dev.off()
