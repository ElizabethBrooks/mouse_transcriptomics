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

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/no_interaction_interaction_nested_random"
dir.create(workingDir)
setwd(workingDir)

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/MackLab/metabolic_adaptation/counted/counts_merged_samples.csv", row.names="gene")

# import grouping factor
targets <- read.csv(file="/Users/bamflappy/MackLab/input_data/experimental_design.csv", row.names="sample")

# set input FDR cutoff
cutFDR <- 0.05

# set LFC cut
#cutFC <- 1.2
cutLFC <- log2(1.2)

# trim the data table to remove summary statistics
countsTable <- head(inputTable, - 5)

# sort columns alphabetically
countsTable <- countsTable[, order(colnames(countsTable))]
targets <- targets[order(rownames(targets)), ]

# check sample naming
setdiff(colnames(countsTable), rownames(targets))
setdiff(rownames(targets), colnames(countsTable))

# filter genes by number of counts
isexpr <- rowSums(cpm(countsTable) > 0.1) >= 5
list <- DGEList(countsTable[isexpr, ])

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


## Dream analysis

# specify parallel processing parameters
param <- SnowParam(4, "SOCK", progressbar = TRUE)

# The variable to be tested must be a fixed effect
form <- ~ 0 + tissue + treatment + sex + mouseline + mouseline:(1|individual) + tissue:mouseline + tissue:sex + treatment:tissue + treatment:mouseline + mouseline:sex

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, targets, BPPARAM = param)

# check the variance partitioning to identify important variables that should be included as fixed or random effects
formVP <- ~ treatment + tissue + sex + mouseline + tissue:mouseline + tissue:sex + treatment:tissue + treatment:mouseline + mouseline:sex
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

# get all results of hypothesis test on treatment
tissueBAT_dge_results <- topTable(fitmm, coef = "tissueBAT", number = nrow(dge))
tissueHYP_dge_results <- topTable(fitmm, coef = "tissueHYP", number = nrow(dge))
tissueLIV_dge_results <- topTable(fitmm, coef = "tissueLIV", number = nrow(dge))
tissueWAT_dge_results <- topTable(fitmm, coef = "tissueWAT", number = nrow(dge))
treatment_dge_results <- topTable(fitmm, coef = "treatmentFood_Restriction", number = nrow(dge))
sex_dge_results <- topTable(fitmm, coef = "sexM", number = nrow(dge))
interactionMANF_dge_results <- topTable(fitmm, coef = "treatmentFood_Restriction:mouselineMANF", number = nrow(dge))
interactionSARA_dge_results <- topTable(fitmm, coef = "treatmentFood_Restriction:mouselineSARA", number = nrow(dge))

# export table of DE genes
tissueBAT_dge_results_tbl <- as_tibble(tissueBAT_dge_results, rownames = "gene")
write.table(tissueBAT_dge_results_tbl, file="tissueBAT.csv", sep=",", row.names=FALSE, quote=FALSE)
tissueHYP_dge_results_tbl <- as_tibble(tissueHYP_dge_results, rownames = "gene")
write.table(tissueHYP_dge_results_tbl, file="tissueHYP.csv", sep=",", row.names=FALSE, quote=FALSE)
tissueLIV_dge_results_tbl <- as_tibble(tissueLIV_dge_results, rownames = "gene")
write.table(tissueLIV_dge_results_tbl, file="tissueLIV.csv", sep=",", row.names=FALSE, quote=FALSE)
tissueWAT_dge_results_tbl <- as_tibble(tissueWAT_dge_results, rownames = "gene")
write.table(tissueWAT_dge_results_tbl, file="tissueWAT.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_dge_results_tbl <- as_tibble(treatment_dge_results, rownames = "gene")
write.table(treatment_dge_results_tbl, file="treatmentFood_Restriction.csv", sep=",", row.names=FALSE, quote=FALSE)
sex_dge_results_tbl <- as_tibble(sex_dge_results, rownames = "gene")
write.table(sex_dge_results_tbl, file="sexM.csv", sep=",", row.names=FALSE, quote=FALSE)
interactionMANF_dge_results_tbl <- as_tibble(interactionMANF_dge_results, rownames = "gene")
write.table(interactionMANF_dge_results_tbl, file="treatmentFood_Restriction_mouselineMANF", sep=",", row.names=FALSE, quote=FALSE)
interactionSARA_dge_results_tbl <- as_tibble(interactionSARA_dge_results, rownames = "gene")
write.table(interactionSARA_dge_results_tbl, file="treatmentFood_Restriction_mouselineSARA", sep=",", row.names=FALSE, quote=FALSE)

# get genes < FDR and LFC cutoffs
tissueBAT_dge_sig <- topTable(fitmm, coef = "tissueBAT", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
tissueHYP_dge_sig <- topTable(fitmm, coef = "tissueHYP", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
tissueLIV_dge_sig <- topTable(fitmm, coef = "tissueLIV", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
tissueWAT_dge_sig <- topTable(fitmm, coef = "tissueWAT", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_dge_sig <- topTable(fitmm, coef = "treatmentFood_Restriction", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
sex_dge_sig <- topTable(fitmm, coef = "sexM", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
interactionMANF_dge_sig <- topTable(fitmm, coef = "treatmentFood_Restriction:mouselineMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
interactionSARA_dge_sig <- topTable(fitmm, coef = "treatmentFood_Restriction:mouselineSARA", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)

# export table of sig DE genes
tissueBAT_dge_sig_tbl <- as_tibble(tissueBAT_dge_sig, rownames = "gene")
tissueBAT_out_file <- paste("tissueBAT", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
tissueBAT_out_file <- paste(tissueBAT_out_file, "csv", sep = ".")
write.table(tissueBAT_dge_sig_tbl, file=tissueBAT_out_file, sep=",", row.names=FALSE, quote=FALSE)
tissueLIV_dge_sig_tbl <- as_tibble(tissueLIV_dge_sig, rownames = "gene")
tissueLIV_out_file <- paste("tissueLIV", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
tissueLIV_out_file <- paste(tissueLIV_out_file, "csv", sep = ".")
write.table(tissueLIV_dge_sig_tbl, file=tissueLIV_out_file, sep=",", row.names=FALSE, quote=FALSE)
tissueHYP_dge_sig_tbl <- as_tibble(tissueHYP_dge_sig, rownames = "gene")
tissueHYP_out_file <- paste("tissueHYP", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
tissueHYP_out_file <- paste(tissueHYP_out_file, "csv", sep = ".")
write.table(tissueBAT_dge_sig_tbl, file=tissueHYP_out_file, sep=",", row.names=FALSE, quote=FALSE)
tissueWAT_dge_sig_tbl <- as_tibble(tissueWAT_dge_sig, rownames = "gene")
tissueWAT_out_file <- paste("tissueWAT", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
tissueWAT_out_file <- paste(tissueWAT_out_file, "csv", sep = ".")
write.table(tissueWAT_dge_sig_tbl, file=tissueWAT_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_dge_sig_tbl <- as_tibble(treatment_dge_sig, rownames = "gene")
treatment_out_file <- paste("treatmentFood_Restriction", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_out_file <- paste(treatment_out_file, "csv", sep = ".")
write.table(treatment_dge_sig_tbl, file=treatment_out_file, sep=",", row.names=FALSE, quote=FALSE)
sex_dge_sig_tbl <- as_tibble(sex_dge_sig, rownames = "gene")
sex_out_file <- paste("sexM", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
sex_out_file <- paste(sex_out_file, "csv", sep = ".")
write.table(sex_dge_sig_tbl, file=sex_out_file, sep=",", row.names=FALSE, quote=FALSE)
interactionMANF_dge_sig_tbl <- as_tibble(interactionMANF_dge_sig, rownames = "gene")
interactionMANF_out_file <- paste("treatmentFood_Restriction_mouselineMANF", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
interactionMANF_out_file <- paste(interactionMANF_out_file, "csv", sep = ".")
write.table(interactionMANF_dge_sig_tbl, file=interactionMANF_out_file, sep=",", row.names=FALSE, quote=FALSE)
interactionSARA_dge_sig_tbl <- as_tibble(interactionSARA_dge_sig, rownames = "gene")
interactionSARA_out_file <- paste("treatmentFood_Restriction_mouselineMANF", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
interactionSARA_out_file <- paste(interactionSARA_out_file, "csv", sep = ".")
write.table(interactionSARA_dge_sig_tbl, file=interactionSARA_out_file, sep=",", row.names=FALSE, quote=FALSE)

# check the number of sig DE genes
nrow(tissueBAT_dge_sig)
nrow(tissueLIV_dge_sig)
nrow(tissueHYP_dge_sig)
nrow(tissueWAT_dge_sig)
nrow(treatment_dge_sig)
nrow(sex_dge_sig)
nrow(interactionMANF_dge_sig)
nrow(interactionSARA_dge_sig)

# add column for identifying direction of DE gene expression
treatment_dge_sig_tbl$colorDE <- plotColors[5]
treatment_dge_sig_tbl$alphaDE <- 0.75
sex_dge_sig_tbl$colorDE <- plotColors[5]
sex_dge_sig_tbl$alphaDE <- 0.75
#interaction_dge_sig_tbl$colorDE <- plotColors[5]
#interaction_dge_sig_tbl$alphaDE <- 0.75

# identify significantly up DE genes
treatment_dge_sig_tbl$colorDE[treatment_dge_sig_tbl$logFC > cutLFC & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
treatment_dge_sig_tbl$alphaDE[treatment_dge_sig_tbl$logFC > cutLFC & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
sex_dge_sig_tbl$colorDE[sex_dge_sig_tbl$logFC > cutLFC & sex_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
sex_dge_sig_tbl$alphaDE[sex_dge_sig_tbl$logFC > cutLFC & sex_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
#interaction_dge_sig_tbl$colorDE[interaction_dge_sig_tbl$logFC > cutLFC & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[4]
#interaction_dge_sig_tbl$alphaDE[interaction_dge_sig_tbl$logFC > cutLFC & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- 1

# identify significantly down DE genes
treatment_dge_sig_tbl$colorDE[treatment_dge_sig_tbl$logFC < (-1*cutLFC) & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
treatment_dge_sig_tbl$alphaDE[treatment_dge_sig_tbl$logFC < (-1*cutLFC) & treatment_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
sex_dge_sig_tbl$colorDE[sex_dge_sig_tbl$logFC < (-1*cutLFC) & sex_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
sex_dge_sig_tbl$alphaDE[sex_dge_sig_tbl$logFC < (-1*cutLFC) & sex_dge_sig_tbl$adj.P.Val < cutFDR] <- 1
#interaction_dge_sig_tbl$colorDE[interaction_dge_sig_tbl$logFC < (-1*cutLFC) & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- plotColors[6]
#interaction_dge_sig_tbl$alphaDE[interaction_dge_sig_tbl$logFC < (-1*cutLFC) & interaction_dge_sig_tbl$adj.P.Val < cutFDR] <- 1

# add column with -log10(adj.P.Val) values
treatment_dge_sig_tbl$negLog10FDR <- -log10(treatment_dge_sig_tbl$adj.P.Val)
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
jpeg("sexM_volcano.jpg")
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
jpeg("sexM_heatmap.jpg")
heatmap(as.matrix(logcountsSubset_sex), margins = c(8, 1), labRow = FALSE)
dev.off()
#jpeg("treatmentUV_groupHT_heatmap.jpg")
#heatmap(as.matrix(logcountsSubset_interaction), margins = c(8, 1), labRow = FALSE)
#dev.off()
