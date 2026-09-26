# inport libraries
library(WGCNA)
library(dplyr)

#The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()

# increase the max size of the vector heap
mem.maxVSize(vsize = "Inf")

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/NAnalysis_22Sep2026/unsigned_power2"
dir.create(workingDir)
setwd(workingDir)

# import grouping factor
factors <- read.csv(file="/Users/bamflappy/MackLab/input_data/experimental_design.csv", row.names="sample")

# retrieve the trait data
traitFile <- "/Users/bamflappy/MackLab/input_data/experimental_design_WGCNA.csv"
allTraits <- read.csv(traitFile, row.names=1)
#allTraits <- allTraits %>% select(-"genotype")

#str(targets)

# read the file
dataFile <- "/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random_tissues/normalizedCounts_logTransformed.csv"
dataInput <- read.csv(file = dataFile, row.names=1)

# transpose each subset
datExpr0 = data = as.data.frame(t(dataInput))
names(datExpr0) = rownames(dataInput)
rownames(datExpr0) = names(dataInput)

#Check the genes across all samples
gsg = goodSamplesGenes(datExpr0, verbose = 3)

# remove the offending genes and samples from the data
if (!gsg$allOK){
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# cluster the samples to see if there are any obvious outliers
sampleTree = hclust(dist(datExpr0), method = "average")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = as.integer(max(sampleTree$height))+1, minSize = 1)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
# filter out samples
datExpr <- datExpr0
datExpr = datExpr0[keepSamples, ]

# retrieve prepared data
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Form a data frame analogous to expression data that will hold the traits
samples = rownames(datExpr)
traitRows = match(samples, rownames(allTraits))
datTraits = allTraits[traitRows,]

# clean up memory
collectGarbage()

# Re-cluster samples
sizeGrWindow(10,7)
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)

# check valid inputs
testCheck <- try(
  plotDendroAndColors(sampleTree2, traitColors),
  silent = TRUE
)
if(class(testCheck) == "try-error"){
  print("FALSE")
}

# Plot the sample dendrogram and the colors underneath.
png(file = "DendroAndColors.png", wi = 9, he = 5, units="in", res=150)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Choose a set of soft-thresholding powers
#powers =  c(seq(from = 1, to = 20, by = 2))
powers = c(c(1:10), seq(from = 12, to=36, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results
cex1 = 0.9
png(file = "SoftPowers.png", wi = 9, he = 5, units="in", res=150)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
abline(h=0.90,col="blue")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = as.integer(nGenes/nSamples)

# set the soft thresholding power
softPower = 2

# determine adjacency
adjacency = adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)

# clean up memory
collectGarbage()

# get distances
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
png(file = "geneClustering.png", wi = 12, he = 9, units="in", res=150)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
png(file = "dynamicTreeCut.png", wi = 8, he = 6, units="in", res=150)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
png(file = "clusteringME.png", wi = 7, he = 6, units="in", res=150)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plot the gene dendrogram again, with the 
# original and merged module colors underneath
png(file = "geneDendro-3.png", wi = 12, he = 9, units="in", res=150)
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkConstruction-stepByStep.RData")

# transpose data
datExpr0 = data = as.data.frame(t(MEs))
names(datExpr0) = rownames(MEs)
rownames(datExpr0) = names(MEs)

# add a column for the row names
datExpr0 <- cbind(gene = rownames(datExpr0), datExpr0)
rownames(datExpr0) <- NULL

# export the expression data as a csv file
write.table(datExpr0, file="eigengeneExpression.csv", sep=",", row.names=FALSE)

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
# create list of module colors mapped to numbers
numMods <- length(unique(moduleColors))
colorTable <- data.frame(
  color = unique(moduleColors),
  number = seq(from = 1, to = numMods, by = 1)
)
# initialize module data frame
resultsTable <- data.frame(
  gene = character(),
  color = character(),
  number = numeric()
)
# match gene IDs with module colors
for(i in 1:numMods){
  gene <- names(datExpr)[moduleColors==colorTable[i,1]]
  color <- rep(colorTable[i,1], length(gene))
  number <- rep(colorTable[i,2], length(gene))
  moduleData <- cbind(gene, color, number)
  resultsTable <- rbind(resultsTable, moduleData)
}

# output table
write.table(resultsTable, "module_genes.csv", sep=",", row.names=FALSE, quote=FALSE)
