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

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/NAnalysis_22Sep2026/unsigned_domesticus"
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

# remove castaneus samples
remove <- rownames(allTraits[grepl("CAST", factors$mouseline),])
dataInput <- select(all_of(dataInput), -c(remove))
allTraits <- allTraits[!rownames(allTraits) %in% remove, ]

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
