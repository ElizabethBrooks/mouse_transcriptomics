#!/usr/bin/env Rscript

# R script to perform DE analysis with DREAM

# turn off scientific notation
options(scipen = 999)

# import libraries
library(variancePartition)
library(edgeR)
library(BiocParallel)
library(tibble)
library(rcartocolor)
library(ggVennDiagram)
library(dplyr)

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/interaction_nested_random"
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

# setup a design matrix
sample_group_pca <- factor(paste(targets$mouseline, targets$treatment, sep="."))

# list sample levels
levels(sample_group_pca)

# setup points and colors for PCA
points <- c(0,0,1,1,2,2)
colors <-  rep(c(plotColors[4], plotColors[6]), 3)

# create a PCA plot with a legend
png("mouseline_treatment_plotPCA.png", units="in", width=6, height=5, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common")
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=3, cex = 0.8)
dev.off()

# setup a design matrix
sample_group_pca <- factor(paste(targets$mouseline, targets$sex, sep="."))

# list sample levels
levels(sample_group_pca)

# setup points and colors for PCA
points <- c(0,0,1,1,2,2)
colors <-  rep(c(plotColors[4], plotColors[6]), 3)

# create a PCA plot with a legend
png("mouseline_sex_plotPCA.png", units="in", width=6, height=5, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common")
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=3)
dev.off()

# setup a design matrix
sample_group_pca <- factor(paste(targets$mouseline, targets$tissue, sep="."))

# list sample levels
levels(sample_group_pca)

# setup points and colors for PCA
points <- c(0,0,0,0,1,1,1,1,2,2,2,2)
#points <- c(15,15,15,15,16,16,16,16,17,17,17,17)
colors <-  rep(c(plotColors[4], plotColors[5], plotColors[6], plotColors[7]), 3)

# create a PCA plot with a legend
png("mouseline_tissue_plotPCA.png", units="in", width=8, height=7, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common")
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=3)
dev.off()

# create a PCA plot with a legend
png("mouseline_tissue_plotPCA_pc3and4.png", units="in", width=8, height=7, res=300)
par(mar=c(4.1, 4.1, 5.1, 0.1), xpd=TRUE)
plotMDS(dge, col=colors[sample_group_pca], pch=points[sample_group_pca], gene.selection="common", dim.plot = c(3, 4))
legend("top", inset=c(0,-0.2), legend=levels(sample_group_pca), pch=points, col=colors, ncol=3)
dev.off()


## Dream analysis

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

# estimate regression coefficients
fitmm <- dream(vobjDream, form, targets)

# apply empirical Bayes shrinkage on linear mixed models
fitmm <- eBayes(fitmm)

# examine design matrix
colnames(fitmm$design)

# treatment: Food_Restriction vs Ad_lib (12-hr food restriction vs. unlimited food, male mice only)
# tissue: LIV (liver), HYP (hypothalamus), BAT (brown adipose), WAT (white adipose)
# mouseline: MANF vs SARA (skinny vs fat)
# sex: M vs F (unlimited food)

# How does food restriction affect gene expression in M. m. domesticus?
## treatmentFood_Restriction
# How do SARA respond differently than MANF to 12-hr food restriction? 
## treatmentFood_Restriction:mouselineMANF
# What genes are expressed differently in M. m. domesticus F compared to M?
## sexF
# What genes are expressed differently in SARA compared to MANF?
## mouselineMANF

# get all results of hypothesis test on treatment
treatment_dge_results <- topTable(fitmm, coef = "treatmentFood_Restriction", number = nrow(dge))
genotype_dge_results <- topTable(fitmm, coef = "mouselineMANF", number = nrow(dge))
sex_dge_results <- topTable(fitmm, coef = "sexF", number = nrow(dge))
interactionMANF_dge_results <- topTable(fitmm, coef = "treatmentFood_Restriction:mouselineMANF", number = nrow(dge))

# export table of DE genes
treatment_dge_results_tbl <- as_tibble(treatment_dge_results, rownames = "gene")
write.table(treatment_dge_results_tbl, file="treatmentFood_Restriction.csv", sep=",", row.names=FALSE, quote=FALSE)
genotype_dge_results_tbl <- as_tibble(genotype_dge_results, rownames = "gene")
write.table(genotype_dge_results_tbl, file="mouselineMANF.csv", sep=",", row.names=FALSE, quote=FALSE)
sex_dge_results_tbl <- as_tibble(sex_dge_results, rownames = "gene")
write.table(sex_dge_results_tbl, file="sexF.csv", sep=",", row.names=FALSE, quote=FALSE)
sex_dge_results_tbl <- as_tibble(sex_dge_results, rownames = "gene")
write.table(sex_dge_results_tbl, file="sexF.csv", sep=",", row.names=FALSE, quote=FALSE)
interactionMANF_dge_results_tbl <- as_tibble(interactionMANF_dge_results, rownames = "gene")
write.table(interactionMANF_dge_results_tbl, file="treatmentFood_Restriction_mouselineMANF", sep=",", row.names=FALSE, quote=FALSE)

# get genes < FDR and LFC cutoffs
treatment_dge_sig <- topTable(fitmm, coef = "treatmentFood_Restriction", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
genotype_dge_sig <- topTable(fitmm, coef = "mouselineMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
sex_dge_sig <- topTable(fitmm, coef = "sexF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
interactionMANF_dge_sig <- topTable(fitmm, coef = "treatmentFood_Restriction:mouselineMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)

# export table of sig DE genes
treatment_dge_sig_tbl <- as_tibble(treatment_dge_sig, rownames = "gene")
treatment_out_file <- paste("treatmentFood_Restriction", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_out_file <- paste(treatment_out_file, "csv", sep = ".")
write.table(treatment_dge_sig_tbl, file=treatment_out_file, sep=",", row.names=FALSE, quote=FALSE)
genotype_dge_sig_tbl <- as_tibble(genotype_dge_sig, rownames = "gene")
genotype_out_file <- paste("mouselineMANF", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
genotype_out_file <- paste(genotype_out_file, "csv", sep = ".")
write.table(genotype_dge_sig_tbl, file=genotype_out_file, sep=",", row.names=FALSE, quote=FALSE)
sex_dge_sig_tbl <- as_tibble(sex_dge_sig, rownames = "gene")
sex_out_file <- paste("sexF", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
sex_out_file <- paste(sex_out_file, "csv", sep = ".")
write.table(sex_dge_sig_tbl, file=sex_out_file, sep=",", row.names=FALSE, quote=FALSE)
interactionMANF_dge_sig_tbl <- as_tibble(interactionMANF_dge_sig, rownames = "gene")
interactionMANF_out_file <- paste("treatmentFood_Restriction_mouselineMANF", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
interactionMANF_out_file <- paste(interactionMANF_out_file, "csv", sep = ".")
write.table(interactionMANF_dge_sig_tbl, file=interactionMANF_out_file, sep=",", row.names=FALSE, quote=FALSE)

# check the number of sig DE genes
nrow(treatment_dge_sig)
nrow(genotype_dge_sig)
nrow(sex_dge_sig)
nrow(interactionMANF_dge_sig)

# add column for identifying direction of DE gene expression
treatment_dge_sig_tbl$colorDE <- plotColors[5]
treatment_dge_sig_tbl$alphaDE <- 0.75
genotype_dge_sig_tbl$colorDE <- plotColors[5]
genotype_dge_sig_tbl$alphaDE <- 0.75
sex_dge_sig_tbl$colorDE <- plotColors[5]
sex_dge_sig_tbl$alphaDE <- 0.75
#interaction_dge_sig_tbl$colorDE <- plotColors[5]
#interaction_dge_sig_tbl$alphaDE <- 0.75

# identify significantly up DE genes
treatment_dge_sig_tbl$colorDE[treatment_dge_sig_tbl$logFC > cutLFC & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
treatment_dge_sig_tbl$alphaDE[treatment_dge_sig_tbl$logFC > cutLFC & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
genotype_dge_sig_tbl$colorDE[genotype_dge_sig_tbl$logFC > cutLFC & genotype_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
genotype_dge_sig_tbl$alphaDE[genotype_dge_sig_tbl$logFC > cutLFC & genotype_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
sex_dge_sig_tbl$colorDE[sex_dge_sig_tbl$logFC > cutLFC & sex_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
sex_dge_sig_tbl$alphaDE[sex_dge_sig_tbl$logFC > cutLFC & sex_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
#interaction_dge_sig_tbl$colorDE[interaction_dge_sig_tbl$logFC > cutLFC & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
#interaction_dge_sig_tbl$alphaDE[interaction_dge_sig_tbl$logFC > cutLFC & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- 1

# identify significantly down DE genes
treatment_dge_sig_tbl$colorDE[treatment_dge_sig_tbl$logFC < (-1*cutLFC) & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
treatment_dge_sig_tbl$alphaDE[treatment_dge_sig_tbl$logFC < (-1*cutLFC) & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
genotype_dge_sig_tbl$colorDE[genotype_dge_sig_tbl$logFC < (-1*cutLFC) & genotype_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
genotype_dge_sig_tbl$alphaDE[genotype_dge_sig_tbl$logFC < (-1*cutLFC) & genotype_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
sex_dge_sig_tbl$colorDE[sex_dge_sig_tbl$logFC < (-1*cutLFC) & sex_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
sex_dge_sig_tbl$alphaDE[sex_dge_sig_tbl$logFC < (-1*cutLFC) & sex_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
#interaction_dge_sig_tbl$colorDE[interaction_dge_sig_tbl$logFC < (-1*cutLFC) & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
#interaction_dge_sig_tbl$alphaDE[interaction_dge_sig_tbl$logFC < (-1*cutLFC) & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- 1

# add column with -log10(adj.P.Val) values
treatment_dge_sig_tbl$negLog10FDR <- -log10(treatment_dge_sig_tbl$adj.P.Val)
genotype_dge_sig_tbl$negLog10FDR <- -log10(genotype_dge_sig_tbl$adj.P.Val)
sex_dge_sig_tbl$negLog10FDR <- -log10(sex_dge_sig_tbl$adj.P.Val)
#interaction_dge_sig_tbl$negLog10FDR <- -log10(interaction_dge_sig_tbl$adj.P.Val)

# create volcano plot
jpeg("treatmentFood_Restriction_volcano.jpg")
ggplot(data=treatment_dge_sig_tbl, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
  geom_point() +
  theme_minimal() +
  scale_color_identity() +
  scale_alpha(guide = 'none') +
  xlab("LFC")
dev.off()
jpeg("mouselineMANF_volcano.jpg")
ggplot(data=genotype_dge_sig_tbl, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
  geom_point() +
  theme_minimal() +
  scale_color_identity() +
  scale_alpha(guide = 'none') +
  xlab("LFC")
dev.off()
jpeg("sexF_volcano.jpg")
ggplot(data=sex_dge_sig_tbl, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
  geom_point() +
  theme_minimal() +
  scale_color_identity() +
  scale_alpha(guide = 'none') +
  xlab("LFC")
dev.off()
#jpeg("treatmentUV_groupHT_volcano.jpg")
#ggplot(data=interaction_dge_sig_tbl, aes(x=logFC, y=negLog10FDR, color = colorDE, alpha = alphaDE)) + 
#  geom_point() +
#  theme_minimal() +
#  scale_color_identity() +
#  scale_alpha(guide = 'none') +
#  xlab("LFC")
#dev.off()

# subset counts table by DE gene set
DGESubset_treatment <- treatment_dge_sig_tbl[!grepl(plotColors[5], treatment_dge_sig_tbl$colorDE),]
DGESubset_treatment.keep <- normListLog$gene %in% DGESubset_treatment$gene
logcountsSubset_treatment <- normListLog[DGESubset_treatment.keep, ]
DGESubset_genotype <- genotype_dge_sig_tbl[!grepl(plotColors[5], genotype_dge_sig_tbl$colorDE),]
DGESubset_genotype.keep <- normListLog$gene %in% DGESubset_genotype$gene
logcountsSubset_genotype <- normListLog[DGESubset_genotype.keep, ]
DGESubset_sex <- sex_dge_sig_tbl[!grepl(plotColors[5], sex_dge_sig_tbl$colorDE),]
DGESubset_sex.keep <- normListLog$gene %in% DGESubset_sex$gene
logcountsSubset_sex <- normListLog[DGESubset_sex.keep, ]
#DGESubset_interaction <- interaction_dge_sig_tbl[!grepl(plotColors[5], interaction_dge_sig_tbl$colorDE),]
#DGESubset_interaction.keep <- normListLog$gene %in% DGESubset_interaction$gene
#logcountsSubset_interaction <- normListLog[DGESubset_interaction.keep, ]

# format for plotting
gene_names <- logcountsSubset_treatment$gene
logcountsSubset_treatment$gene <- NULL
rownames(logcountsSubset_treatment) <- gene_names
gene_names <- logcountsSubset_genotype$gene
logcountsSubset_genotype$gene <- NULL
rownames(logcountsSubset_genotype) <- gene_names
gene_names <- logcountsSubset_sex$gene
logcountsSubset_sex$gene <- NULL
rownames(logcountsSubset_sex) <- gene_names
#gene_names <- logcountsSubset_interaction$gene
#logcountsSubset_interaction$gene <- NULL
#rownames(logcountsSubset_interaction) <- gene_names

# heatmap of results
jpeg("treatmentFood_Restriction_heatmap.jpg")
heatmap(as.matrix(logcountsSubset_treatment), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("mouselineMANF_heatmap.jpg")
heatmap(as.matrix(logcountsSubset_genotype), margins = c(8, 1), labRow = FALSE)
dev.off()
jpeg("sexF_heatmap.jpg")
heatmap(as.matrix(logcountsSubset_sex), margins = c(8, 1), labRow = FALSE)
dev.off()
#jpeg("treatmentUV_groupHT_heatmap.jpg")
#heatmap(as.matrix(logcountsSubset_interaction), margins = c(8, 1), labRow = FALSE)
#dev.off()
