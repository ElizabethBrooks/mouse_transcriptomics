#!/usr/bin/env Rscript

# R script to perform DE analysis with DREAM

# turn off scientific notation
options(scipen = 999)

# import libraries
library(edgeR)
library(tibble)
library(rcartocolor)
library(dplyr)
library(scales)

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])
plotColorSubset <- c(plotColors[3], plotColors[4], plotColors[5], plotColors[6])
show_col(plotColors)
show_col(plotColorSubset)

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/data_exploration_domesticus"
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

# remove castaneus samples
remove <- rownames(factors[grepl("CAST", factors$mouseline),])
countsTable <- select(all_of(countsTable), -c(remove))
factors <- factors[!grepl("CAST", factors$mouseline),]

# convert the grouping data into factors 
targets <- as.data.frame(lapply(factors, as.factor))
rownames(targets) <- rownames(factors)

# specify the reference level
targets$treatment <- relevel(targets$treatment, ref = "Ad_lib")
targets$mouseline <- relevel(targets$mouseline, ref = "SARA")
targets$tissue <- relevel(targets$tissue, ref = "BAT")
targets$sex <- relevel(targets$sex, ref = "Male")

# sort samples alphabetically
countsTable <- countsTable[, order(colnames(countsTable))]
targets <- targets[order(rownames(targets)), ]

# check sample naming
setdiff(colnames(countsTable), rownames(targets))
setdiff(rownames(targets), colnames(countsTable))

# there are at least 2 samples per individual
smallestGroupSize <- 2

# keep only rows that have a count of at least 10 for a minimal number of samples
keep <- rowSums(cpm(countsTable) >= 10) >= smallestGroupSize
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

# To-do: fix PCA plots

# setup a design matrix
sample_group_pca <- factor(paste(targets$mouseline, targets$treatment, sep="."))

# list sample levels
levels(sample_group_pca)

# setup points and colors for PCA
points <- c(0,0,1,1)
colors <-  rep(c(plotColors[4], plotColors[6]), 2)

# create a PCA plot with a legend
png("mouseline_treatment_plotPCA.png", units="in", width=6, height=5, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common")
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=2, cex = 0.8)
dev.off()

# create a PCA plot with a legend
png("mouseline_treatment_plotPCA_pc3and4.png", units="in", width=6, height=5, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common", dim.plot = c(3, 4))
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=2, cex = 0.8)
dev.off()

# setup a design matrix
sample_group_pca <- factor(paste(targets$mouseline, targets$sex, sep="."))

# list sample levels
levels(sample_group_pca)

# setup points and colors for PCA
points <- c(0,0,1,1)
colors <-  rep(c(plotColors[4], plotColors[6]), 2)

# create a PCA plot with a legend
png("mouseline_sex_plotPCA.png", units="in", width=6, height=5, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common")
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=2)
dev.off()

# create a PCA plot with a legend
png("mouseline_sex_plotPCA_pc3and4.png", units="in", width=6, height=5, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common", dim.plot = c(3, 4))
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=2)
dev.off()

# setup a design matrix
sample_group_pca <- factor(paste(targets$mouseline, targets$tissue, sep="."))

# list sample levels
levels(sample_group_pca)

# setup points and colors for PCA
points <- c(0,0,0,0,1,1,1,1)
#points <- c(15,15,15,15,16,16,16,16,17,17,17,17)
colors <-  rep(c(plotColorSubset), 2)

# create a PCA plot with a legend
png("mouseline_tissue_plotPCA.png", units="in", width=8, height=7, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common")
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=2)
dev.off()

# create a PCA plot with a legend
png("mouseline_tissue_plotPCA_pc3and4.png", units="in", width=8, height=7, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common", dim.plot = c(3, 4))
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=2)
dev.off()
