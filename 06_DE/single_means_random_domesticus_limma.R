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

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random_domesticus"
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
  SARAvsMANF=
    eval(parse(text=paste(factor_SARA, collapse = "+")))/length(factor_SARA) -
    eval(parse(text=paste(factor_MANF, collapse = "+")))/length(factor_MANF), 
  M_SARAvsMANF=
    eval(parse(text=paste(factor_SARA_M, collapse = "+")))/length(factor_SARA_M) -
    eval(parse(text=paste(factor_MANF_M, collapse = "+")))/length(factor_MANF_M), 
  M_Ad_SARAvsMANF=
    eval(parse(text=paste(factor_SARA_M_Ad, collapse = "+")))/length(factor_SARA_M_Ad) -
    eval(parse(text=paste(factor_MANF_M_Ad, collapse = "+")))/length(factor_MANF_M_Ad), 
  AdvsFood=
    eval(parse(text=paste(factor_Ad, collapse = "+")))/length(factor_Ad) -
    eval(parse(text=paste(factor_Food, collapse = "+")))/length(factor_Food),
  MvsF=
    eval(parse(text=paste(factor_M, collapse = "+")))/length(factor_M) -
    eval(parse(text=paste(factor_F, collapse = "+")))/length(factor_F),
  Ad_MvsF=
    eval(parse(text=paste(factor_M_Ad, collapse = "+")))/length(factor_M_Ad) -
    eval(parse(text=paste(factor_F, collapse = "+")))/length(factor_F),
  SARA_MvsF=
    eval(parse(text=paste(factor_SARA_M, collapse = "+")))/length(factor_SARA_M) -
    eval(parse(text=paste(factor_SARA_F, collapse = "+")))/length(factor_SARA_F),
  MANF_MvsF=
    eval(parse(text=paste(factor_MANF_M, collapse = "+")))/length(factor_MANF_M) -
    eval(parse(text=paste(factor_MANF_F, collapse = "+")))/length(factor_MANF_F),
  SARA_Ad_MvsF=
    eval(parse(text=paste(factor_SARA_M_Ad, collapse = "+")))/length(factor_SARA_M_Ad) -
    eval(parse(text=paste(factor_SARA_F, collapse = "+")))/length(factor_SARA_F),
  MANF_Ad_MvsF=
    eval(parse(text=paste(factor_MANF_M_Ad, collapse = "+")))/length(factor_MANF_M_Ad) -
    eval(parse(text=paste(factor_MANF_F, collapse = "+")))/length(factor_MANF_F),
  SARAvsMANF_MvsF=
    (eval(parse(text=paste(factor_SARA_M, collapse = "+")))/length(factor_SARA_M) -
    eval(parse(text=paste(factor_SARA_F, collapse = "+")))/length(factor_SARA_F)) -
    (eval(parse(text=paste(factor_MANF_M, collapse = "+")))/length(factor_MANF_M) -
       eval(parse(text=paste(factor_MANF_F, collapse = "+")))/length(factor_MANF_F)),
  Ad_SARAvsMANF_MvsF=
    (eval(parse(text=paste(factor_SARA_M_Ad, collapse = "+")))/length(factor_SARA_M_Ad) -
       eval(parse(text=paste(factor_SARA_F, collapse = "+")))/length(factor_SARA_F)) -
    (eval(parse(text=paste(factor_MANF_M_Ad, collapse = "+")))/length(factor_MANF_M_Ad) -
       eval(parse(text=paste(factor_MANF_F, collapse = "+")))/length(factor_MANF_F)),
  M_AdvsFood=
    eval(parse(text=paste(factor_M_Ad, collapse = "+")))/length(factor_M_Ad) -
    eval(parse(text=paste(factor_M_Food, collapse = "+")))/length(factor_M_Food),
  SARA_AdvsFood=
    eval(parse(text=paste(factor_SARA_Ad, collapse = "+")))/length(factor_SARA_Ad) -
    eval(parse(text=paste(factor_SARA_Food, collapse = "+")))/length(factor_SARA_Food),
  MANF_AdvsFood=
    eval(parse(text=paste(factor_MANF_Ad, collapse = "+")))/length(factor_MANF_Ad) -
    eval(parse(text=paste(factor_MANF_Food, collapse = "+")))/length(factor_MANF_Food),
  SARA_M_AdvsFood=
    eval(parse(text=paste(factor_SARA_M_Ad, collapse = "+")))/length(factor_SARA_M_Ad) -
    eval(parse(text=paste(factor_SARA_M_Food, collapse = "+")))/length(factor_SARA_M_Food),
  MANF_M_AdvsFood=
    eval(parse(text=paste(factor_MANF_M_Ad, collapse = "+")))/length(factor_MANF_M_Ad) -
    eval(parse(text=paste(factor_MANF_M_Food, collapse = "+")))/length(factor_MANF_M_Food),
  SARAvsMANF_AdvsFood=
    (eval(parse(text=paste(factor_SARA_Ad, collapse = "+")))/length(factor_SARA_Ad) -
       eval(parse(text=paste(factor_SARA_Food, collapse = "+")))/length(factor_SARA_Food)) - 
    (eval(parse(text=paste(factor_MANF_Ad, collapse = "+")))/length(factor_MANF_Ad) -
       eval(parse(text=paste(factor_MANF_Food, collapse = "+")))/length(factor_MANF_Food)),
  M_SARAvsMANF_AdvsFood=
    (eval(parse(text=paste(factor_SARA_M_Ad, collapse = "+")))/length(factor_SARA_M_Ad) -
       eval(parse(text=paste(factor_SARA_M_Food, collapse = "+")))/length(factor_SARA_M_Food)) - 
    (eval(parse(text=paste(factor_MANF_M_Ad, collapse = "+")))/length(factor_MANF_M_Ad) -
       eval(parse(text=paste(factor_MANF_M_Food, collapse = "+")))/length(factor_MANF_M_Food)),
  BAT_SARAvsMANF=
    eval(parse(text=paste(factor_SARA_BAT, collapse = "+")))/length(factor_SARA_BAT) -
    eval(parse(text=paste(factor_MANF_BAT, collapse = "+")))/length(factor_MANF_BAT),
  WAT_SARAvsMANF=
    eval(parse(text=paste(factor_SARA_WAT, collapse = "+")))/length(factor_SARA_WAT) -
    eval(parse(text=paste(factor_MANF_WAT, collapse = "+")))/length(factor_MANF_WAT),
  HYP_SARAvsMANF=
    eval(parse(text=paste(factor_SARA_HYP, collapse = "+")))/length(factor_SARA_HYP) -
    eval(parse(text=paste(factor_MANF_HYP, collapse = "+")))/length(factor_MANF_HYP),
  LIV_SARAvsMANF=
    eval(parse(text=paste(factor_SARA_LIV, collapse = "+")))/length(factor_SARA_LIV) -
    eval(parse(text=paste(factor_MANF_LIV, collapse = "+")))/length(factor_MANF_LIV),
  SARA_BAT_AdvsFood=
    eval(parse(text=paste(factor_SARA_BAT_Ad, collapse = "+")))/length(factor_SARA_BAT_Ad) -
    eval(parse(text=paste(factor_SARA_BAT_Food, collapse = "+")))/length(factor_SARA_BAT_Food),
  SARA_WAT_AdvsFood=
    eval(parse(text=paste(factor_SARA_WAT_Ad, collapse = "+")))/length(factor_SARA_WAT_Ad) -
    eval(parse(text=paste(factor_SARA_WAT_Food, collapse = "+")))/length(factor_SARA_WAT_Food),
  SARA_HYP_AdvsFood=
    eval(parse(text=paste(factor_SARA_HYP_Ad, collapse = "+")))/length(factor_SARA_HYP_Ad) -
    eval(parse(text=paste(factor_SARA_HYP_Food, collapse = "+")))/length(factor_SARA_HYP_Food),
  SARA_LIV_AdvsFood=
    eval(parse(text=paste(factor_SARA_LIV_Ad, collapse = "+")))/length(factor_SARA_LIV_Ad) -
    eval(parse(text=paste(factor_SARA_LIV_Food, collapse = "+")))/length(factor_SARA_LIV_Food),
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

# What genes are expressed differently in SARA compared to MANF?
## SARAvsMANF
# What genes are expressed differently in male SARA compared to MANF?
## M_SARAvsMANF
# What genes are expressed differently in male SARA compared to MANF when not food restricted?
## M_Ad_SARAvsMANF
# How does 12-hr food restriction affect gene expression in M. m.?
## AdvsFood
# How does 12-hr food restriction affect gene expression in M. m. M?
## M_AdvsFood
# How do SARA respond to 12-hr food restriction? 
## SARA_AdvsFood
# How do SARA M respond to 12-hr food restriction? 
## SARA_M_AdvsFood
# How do MANF respond to 12-hr food restriction? 
## MANF_AdvsFood
# How do MANF M respond to 12-hr food restriction? 
## MANF_M_AdvsFood
# How do SARA respond differently than MANF to 12-hr food restriction? 
## SARAvsMANF_AdvsFood
# How do SARA M respond differently than MANF M to 12-hr food restriction? 
## M_SARAvsMANF_AdvsFood
# How do SARA tissues respond to 12-hr food restriction?
## SARA_BAT_AdvsFood
## SARA_WAT_AdvsFood
## SARA_HYP_AdvsFood
## SARA_LIV_AdvsFood
# How do MANF tissues respond to 12-hr food restriction?
## MANF_BAT_AdvsFood
## MANF_WAT_AdvsFood
## MANF_HYP_AdvsFood
## MANF_LIV_AdvsFood
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
# What genes are expressed differently in SARA tissues compared to MANF tissues?
## BAT_SARAvsMANF
## WAT_SARAvsMANF
## HYP_SARAvsMANF
## LIV_SARAvsMANF
# What genes are expressed differently in M. m. F compared to M?
## MvsF
# What genes are expressed differently in SARA F compared to M?
## SARA_MvsF
# What genes are expressed differently in MANF F compared to M?
## MANF_MvsF
# What genes are expressed differently in M. m. F compared to M when not food restricted?
## Ad_MvsF
# What genes are expressed differently in SARA F compared to M when not food restricted?
## SARA_Ad_MvsF
# What genes are expressed differently in MANF F compared to M when not food restricted?
## MANF_Ad_MvsF
# What genes are expressed differently in SARA vs MANF F compared to M?
## SARAvsMANF_MvsF
# What genes are expressed differently in SARA vs MANF F compared to M when not food restricted?
## Ad_SARAvsMANF_MvsF
# How do SARA tissues respond differently than MANF tissues to 12-hr food restriction?
# What genes are expressed differently in M. m. domesticus F tissues compared to M tissues?
# What genes are expressed differently in SARA F tissues compared to M tissues?
# What genes are expressed differently in MANF F tissues compared to M tissues?

# get genes < FDR and LFC cutoffs
# genotype
genotype_dge_sig <- topTable(fitCont, coef = "SARAvsMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
genotype_M_dge_sig <- topTable(fitCont, coef = "M_SARAvsMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
genotype_M_Ad_dge_sig <- topTable(fitCont, coef = "M_Ad_SARAvsMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
# treatment overall and among species
treatment_dge_sig <- topTable(fitCont, coef = "AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_M_dge_sig <- topTable(fitCont, coef = "M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_dge_sig <- topTable(fitCont, coef = "SARA_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_dge_sig <- topTable(fitCont, coef = "MANF_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_M_dge_sig <- topTable(fitCont, coef = "SARA_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_M_dge_sig <- topTable(fitCont, coef = "MANF_M_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_dge_sig <- topTable(fitCont, coef = "SARAvsMANF_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_M_dge_sig <- topTable(fitCont, coef = "M_SARAvsMANF_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
# treatment among tissues
treatment_SARA_BAT_dge_sig <- topTable(fitCont, coef = "SARA_BAT_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_WAT_dge_sig <- topTable(fitCont, coef = "SARA_WAT_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_HYP_dge_sig <- topTable(fitCont, coef = "SARA_HYP_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_LIV_dge_sig <- topTable(fitCont, coef = "SARA_LIV_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_BAT_dge_sig <- topTable(fitCont, coef = "MANF_BAT_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_WAT_dge_sig <- topTable(fitCont, coef = "MANF_WAT_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_HYP_dge_sig <- topTable(fitCont, coef = "MANF_HYP_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_MANF_LIV_dge_sig <- topTable(fitCont, coef = "MANF_LIV_AdvsFood", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_BAT_dge_sig <- topTable(fitCont, coef = "BAT_SARAvsMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_WAT_dge_sig <- topTable(fitCont, coef = "WAT_SARAvsMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_HYP_dge_sig <- topTable(fitCont, coef = "HYP_SARAvsMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
treatment_SARA_MANF_LIV_dge_sig <- topTable(fitCont, coef = "LIV_SARAvsMANF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
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
# sex overall and among species
sex_dge_sig <- topTable(fitCont, coef = "MvsF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
sex_Ad_dge_sig <- topTable(fitCont, coef = "Ad_MvsF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
sex_SARA_dge_sig <- topTable(fitCont, coef = "SARA_MvsF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
sex_MANF_dge_sig <- topTable(fitCont, coef = "MANF_MvsF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
sex_SARA_Ad_dge_sig <- topTable(fitCont, coef = "SARA_Ad_MvsF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
sex_MANF_Ad_dge_sig <- topTable(fitCont, coef = "MANF_Ad_MvsF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
sex_SARA_MANF_dge_sig <- topTable(fitCont, coef = "SARAvsMANF_MvsF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)
sex_SARA_MANF_Ad_dge_sig <- topTable(fitCont, coef = "Ad_SARAvsMANF_MvsF", number = nrow(dge), p.value = cutFDR, lfc = cutLFC)

# check the number of sig DE genes
# genotype
#nrow(genotype_dge_sig)
#nrow(genotype_M_dge_sig)
nrow(genotype_M_Ad_dge_sig)
# treatment overall and among species
#nrow(treatment_dge_sig)
nrow(treatment_M_dge_sig)
#nrow(treatment_SARA_dge_sig)
#nrow(treatment_MANF_dge_sig)
nrow(treatment_SARA_M_dge_sig)
nrow(treatment_MANF_M_dge_sig)
#nrow(treatment_SARA_MANF_dge_sig)
nrow(treatment_SARA_MANF_M_dge_sig)
# treatment among tissues
#nrow(treatment_SARA_BAT_dge_sig)
#nrow(treatment_SARA_WAT_dge_sig)
#nrow(treatment_SARA_HYP_dge_sig)
#nrow(treatment_SARA_LIV_dge_sig)
#nrow(treatment_MANF_BAT_dge_sig)
#nrow(treatment_MANF_WAT_dge_sig)
#nrow(treatment_MANF_HYP_dge_sig)
#nrow(treatment_MANF_LIV_dge_sig)
#nrow(treatment_SARA_MANF_BAT_dge_sig)
#nrow(treatment_SARA_MANF_WAT_dge_sig)
#nrow(treatment_SARA_MANF_HYP_dge_sig)
#nrow(treatment_SARA_MANF_LIV_dge_sig)
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
# sex overall and among species
#nrow(sex_dge_sig)
#nrow(sex_Ad_dge_sig)
#nrow(sex_SARA_dge_sig)
#nrow(sex_MANF_dge_sig)
nrow(sex_SARA_Ad_dge_sig)
nrow(sex_MANF_Ad_dge_sig)
#nrow(sex_SARA_MANF_dge_sig)
nrow(sex_SARA_MANF_Ad_dge_sig)

# export tables of sig DE genes
# genotype: genotype_M_Ad_dge_sig
genotype_M_Ad_dge_sig_tbl <- as_tibble(genotype_M_Ad_dge_sig, rownames = "gene")
genotype_M_Ad_out_file <- paste("genotype_M_Ad_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
genotype_M_Ad_out_file <- paste(genotype_M_Ad_out_file, "csv", sep = ".")
write.table(genotype_M_Ad_dge_sig_tbl, file=genotype_M_Ad_out_file, sep=",", row.names=FALSE, quote=FALSE)
# treatment: treatment_M_dge_sig, treatment_SARA_M_dge_sig, treatment_MANF_M_dge_sig, treatment_SARA_MANF_M_dge_sig
treatment_M_dge_sig_tbl <- as_tibble(treatment_M_dge_sig, rownames = "gene")
treatment_M_out_file <- paste("treatment_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_M_out_file <- paste(treatment_M_out_file, "csv", sep = ".")
write.table(treatment_M_dge_sig_tbl, file=treatment_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_M_dge_sig_tbl <- as_tibble(treatment_SARA_M_dge_sig, rownames = "gene")
treatment_SARA_M_out_file <- paste("treatment_SARA_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_M_out_file <- paste(treatment_SARA_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_M_dge_sig_tbl, file=treatment_SARA_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_M_dge_sig_tbl <- as_tibble(treatment_MANF_M_dge_sig, rownames = "gene")
treatment_MANF_M_out_file <- paste("treatment_MANF_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_MANF_M_out_file <- paste(treatment_MANF_M_out_file, "csv", sep = ".")
write.table(treatment_MANF_M_dge_sig_tbl, file=treatment_MANF_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_M_dge_sig_tbl <- as_tibble(treatment_SARA_MANF_M_dge_sig, rownames = "gene")
treatment_SARA_MANF_M_out_file <- paste("treatment_SARA_MANF_M_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
treatment_SARA_MANF_M_out_file <- paste(treatment_SARA_MANF_M_out_file, "csv", sep = ".")
write.table(treatment_SARA_MANF_M_dge_sig_tbl, file=treatment_SARA_MANF_M_out_file, sep=",", row.names=FALSE, quote=FALSE)
# tissue SARA: treatment_SARA_BAT_M_dge_sig, treatment_SARA_WAT_M_dge_sig, treatment_SARA_HYP_M_dge_sig, treatment_SARA_LIV_M_dge_sig
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
# tissue MANF: treatment_MANF_BAT_M_dge_sig, treatment_MANF_WAT_M_dge_sig, treatment_MANF_HYP_M_dge_sig, treatment_MANF_LIV_M_dge_sig
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
# tissue SARA vs MANF: treatment_SARA_MANF_BAT_M_dge_sig, treatment_SARA_MANF_WAT_M_dge_sig, treatment_SARA_MANF_HYP_M_dge_sig, treatment_SARA_MANF_LIV_M_dge_sig
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
# sex: sex_SARA_Ad_dge_sig, sex_MANF_Ad_dge_sig, sex_SARA_MANF_Ad_dge_sig
sex_SARA_Ad_dge_sig_tbl <- as_tibble(sex_SARA_Ad_dge_sig, rownames = "gene")
sex_SARA_Ad_out_file <- paste("sex_SARA_Ad_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
sex_SARA_Ad_out_file <- paste(sex_SARA_Ad_out_file, "csv", sep = ".")
write.table(sex_SARA_Ad_dge_sig_tbl, file=sex_SARA_Ad_out_file, sep=",", row.names=FALSE, quote=FALSE)
sex_MANF_Ad_dge_sig_tbl <- as_tibble(sex_MANF_Ad_dge_sig, rownames = "gene")
sex_MANF_Ad_out_file <- paste("sex_MANF_Ad_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
sex_MANF_Ad_out_file <- paste(sex_MANF_Ad_out_file, "csv", sep = ".")
write.table(sex_MANF_Ad_dge_sig_tbl, file=sex_MANF_Ad_out_file, sep=",", row.names=FALSE, quote=FALSE)
sex_SARA_MANF_Ad_dge_sig_tbl <- as_tibble(sex_SARA_MANF_Ad_dge_sig, rownames = "gene")
sex_SARA_MANF_Ad_out_file <- paste("sex_SARA_MANF_Ad_dge_sig", "FDR", cutFDR, "LFC", cutLFC, sep = "_")
sex_SARA_MANF_Ad_out_file <- paste(sex_SARA_MANF_Ad_out_file, "csv", sep = ".")
write.table(sex_SARA_MANF_Ad_dge_sig_tbl, file=sex_SARA_MANF_Ad_out_file, sep=",", row.names=FALSE, quote=FALSE)

# get all results of hypothesis tests
# genotype: genotype_M_Ad_dge_sig
genotype_M_Ad_dge_results <- topTable(fitCont, coef = "M_Ad_SARAvsMANF", number = nrow(dge))
# treatment: treatment_M_dge_sig, treatment_SARA_M_dge_sig, treatment_MANF_M_dge_sig, treatment_SARA_MANF_M_dge_sig
treatment_M_dge_results <- topTable(fitCont, coef = "M_AdvsFood", number = nrow(dge))
treatment_SARA_M_dge_results <- topTable(fitCont, coef = "SARA_M_AdvsFood", number = nrow(dge))
treatment_MANF_M_dge_results <- topTable(fitCont, coef = "MANF_M_AdvsFood", number = nrow(dge))
treatment_SARA_MANF_M_dge_results <- topTable(fitCont, coef = "M_SARAvsMANF_AdvsFood", number = nrow(dge))
# tissue SARA: treatment_SARA_BAT_M_dge_sig, treatment_SARA_WAT_M_dge_sig, treatment_SARA_HYP_M_dge_sig, treatment_SARA_LIV_M_dge_sig
treatment_SARA_BAT_M_dge_results <- topTable(fitCont, coef = "SARA_BAT_M_AdvsFood", number = nrow(dge))
treatment_SARA_WAT_M_dge_results <- topTable(fitCont, coef = "SARA_WAT_M_AdvsFood", number = nrow(dge))
treatment_SARA_HYP_M_dge_results <- topTable(fitCont, coef = "SARA_HYP_M_AdvsFood", number = nrow(dge))
treatment_SARA_LIV_M_dge_results <- topTable(fitCont, coef = "SARA_LIV_M_AdvsFood", number = nrow(dge))
# tissue MANF: treatment_MANF_BAT_M_dge_sig, treatment_MANF_WAT_M_dge_sig, treatment_MANF_HYP_M_dge_sig, treatment_MANF_LIV_M_dge_sig
treatment_MANF_BAT_M_dge_results <- topTable(fitCont, coef = "SARA_BAT_M_AdvsFood", number = nrow(dge))
treatment_MANF_WAT_M_dge_results <- topTable(fitCont, coef = "SARA_WAT_M_AdvsFood", number = nrow(dge))
treatment_MANF_HYP_M_dge_results <- topTable(fitCont, coef = "SARA_HYP_M_AdvsFood", number = nrow(dge))
treatment_MANF_LIV_M_dge_results <- topTable(fitCont, coef = "SARA_LIV_M_AdvsFood", number = nrow(dge))
# tissue SARA vs MANF: treatment_SARA_MANF_BAT_M_dge_sig, treatment_SARA_MANF_WAT_M_dge_sig, treatment_SARA_MANF_HYP_M_dge_sig, treatment_SARA_MANF_LIV_M_dge_sig
treatment_SARA_MANF_BAT_M_dge_results <- topTable(fitCont, coef = "SARAvsMANF_BAT_M_AdvsFood", number = nrow(dge))
treatment_SARA_MANF_WAT_M_dge_results <- topTable(fitCont, coef = "SARAvsMANF_WAT_M_AdvsFood", number = nrow(dge))
treatment_SARA_MANF_HYP_M_dge_results <- topTable(fitCont, coef = "SARAvsMANF_HYP_M_AdvsFood", number = nrow(dge))
treatment_SARA_MANF_LIV_M_dge_results <- topTable(fitCont, coef = "SARAvsMANF_LIV_M_AdvsFood", number = nrow(dge))
# sex: sex_SARA_Ad_dge_sig, sex_MANF_Ad_dge_sig, sex_SARA_MANF_Ad_dge_sig
sex_SARA_Ad_dge_results <- topTable(fitCont, coef = "SARA_Ad_MvsF", number = nrow(dge))
sex_MANF_Ad_dge_results <- topTable(fitCont, coef = "MANF_Ad_MvsF", number = nrow(dge))
sex_SARA_MANF_Ad_dge_results <- topTable(fitCont, coef = "Ad_SARAvsMANF_MvsF", number = nrow(dge))

# export table of DE genes
# genotype: genotype_M_Ad_dge_sig
genotype_M_Ad_dge_results_tbl <- as_tibble(genotype_M_Ad_dge_results, rownames = "gene")
write.table(genotype_M_Ad_dge_results_tbl, file="genotype_M_Ad_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
# treatment: treatment_M_dge_sig, treatment_SARA_M_dge_sig, treatment_MANF_M_dge_sig, treatment_SARA_MANF_M_dge_sig
treatment_M_dge_results_tbl <- as_tibble(treatment_M_dge_results, rownames = "gene")
write.table(treatment_M_dge_results_tbl, file="treatment_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_M_dge_results_tbl <- as_tibble(treatment_SARA_M_dge_results, rownames = "gene")
write.table(treatment_SARA_M_dge_results_tbl, file="treatment_SARA_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_M_dge_results_tbl <- as_tibble(treatment_MANF_M_dge_results, rownames = "gene")
write.table(treatment_MANF_M_dge_results_tbl, file="treatment_MANF_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_M_dge_results_tbl <- as_tibble(treatment_SARA_MANF_M_dge_results, rownames = "gene")
write.table(treatment_SARA_MANF_M_dge_results_tbl, file="treatment_SARA_MANF_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
# tissue SARA: treatment_SARA_BAT_M_dge_sig, treatment_SARA_WAT_M_dge_sig, , treatment_SARA_LIV_M_dge_sig
treatment_SARA_BAT_M_dge_results_tbl <- as_tibble(treatment_SARA_BAT_M_dge_results, rownames = "gene")
write.table(treatment_SARA_BAT_M_dge_results_tbl, file="treatment_SARA_BAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_WAT_M_dge_results_tbl <- as_tibble(treatment_SARA_WAT_M_dge_results, rownames = "gene")
write.table(treatment_SARA_WAT_M_dge_results_tbl, file="treatment_SARA_WAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_HYP_M_dge_results_tbl <- as_tibble(treatment_SARA_HYP_M_dge_results, rownames = "gene")
write.table(treatment_SARA_HYP_M_dge_results_tbl, file="treatment_SARA_HYP_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_LIV_M_dge_results_tbl <- as_tibble(treatment_SARA_LIV_M_dge_results, rownames = "gene")
write.table(treatment_SARA_LIV_M_dge_results_tbl, file="treatment_SARA_LIV_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
# tissue MANF: treatment_MANF_BAT_M_dge_sig, treatment_MANF_WAT_M_dge_sig, , treatment_MANF_LIV_M_dge_sig
treatment_MANF_BAT_M_dge_results_tbl <- as_tibble(treatment_MANF_BAT_M_dge_results, rownames = "gene")
write.table(treatment_MANF_BAT_M_dge_results_tbl, file="treatment_MANF_BAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_WAT_M_dge_results_tbl <- as_tibble(treatment_MANF_WAT_M_dge_results, rownames = "gene")
write.table(treatment_MANF_WAT_M_dge_results_tbl, file="treatment_MANF_WAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_HYP_M_dge_results_tbl <- as_tibble(treatment_MANF_HYP_M_dge_results, rownames = "gene")
write.table(treatment_MANF_HYP_M_dge_results_tbl, file="treatment_MANF_HYP_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_MANF_LIV_M_dge_results_tbl <- as_tibble(treatment_MANF_LIV_M_dge_results, rownames = "gene")
write.table(treatment_MANF_LIV_M_dge_results_tbl, file="treatment_MANF_LIV_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
# tissue SARA vs MANF: treatment_SARA_MANF_BAT_M_dge_sig, treatment_SARA_MANF_WAT_M_dge_sig, , treatment_SARA_MANF_LIV_M_dge_sig
treatment_SARA_MANF_BAT_M_dge_results_tbl <- as_tibble(treatment_SARA_MANF_BAT_M_dge_results, rownames = "gene")
write.table(treatment_SARA_MANF_BAT_M_dge_results_tbl, file="treatment_SARA_MANF_BAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_WAT_M_dge_results_tbl <- as_tibble(treatment_SARA_MANF_WAT_M_dge_results, rownames = "gene")
write.table(treatment_SARA_MANF_WAT_M_dge_results_tbl, file="treatment_SARA_MANF_WAT_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_HYP_M_dge_results_tbl <- as_tibble(treatment_SARA_MANF_HYP_M_dge_results, rownames = "gene")
write.table(treatment_SARA_MANF_HYP_M_dge_results_tbl, file="treatment_SARA_MANF_HYP_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
treatment_SARA_MANF_LIV_M_dge_results_tbl <- as_tibble(treatment_SARA_MANF_LIV_M_dge_results, rownames = "gene")
write.table(treatment_SARA_MANF_LIV_M_dge_results_tbl, file="treatment_SARA_MANF_LIV_M_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
# sex: sex_SARA_Ad_dge_sig, sex_MANF_Ad_dge_sig, sex_SARA_MANF_Ad_dge_sig
sex_SARA_Ad_dge_results_tbl <- as_tibble(sex_SARA_Ad_dge_results, rownames = "gene")
write.table(sex_SARA_Ad_dge_results_tbl, file="sex_SARA_Ad_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
sex_MANF_Ad_dge_results_tbl <- as_tibble(sex_MANF_Ad_dge_results, rownames = "gene")
write.table(sex_MANF_Ad_dge_results_tbl, file="sex_MANF_Ad_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
sex_SARA_MANF_Ad_dge_results_tbl <- as_tibble(sex_SARA_MANF_Ad_dge_results, rownames = "gene")
write.table(sex_SARA_MANF_Ad_dge_results_tbl, file="sex_SARA_MANF_Ad_dge_results.csv", sep=",", row.names=FALSE, quote=FALSE)
