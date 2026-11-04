#!/usr/bin/env Rscript

# script to plot bulk ATAC fragment size distributions
# usage: RScript fragment_sizes_ATAC.R bam_file sample_label

# install the package
#if (!require("BiocManager", quietly = TRUE)) {
#    install.packages("BiocManager")
#}
#BiocManager::install("ATACseqQC")

# load the package
library(ATACseqQC)

# fetch trailing command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/ATAC_analysis_21Oct2026/evaluation"
dir.create(workingDir)
setwd(workingDir)

# path to the input BAM file
#bam_file <- "/Users/bamflappy/MackLab/metabolic_adaptation/sorted_coordinate_ATAC/25164FL-06-01-01_S347/accepted_hits.bam"
bam_file <- args[1]

# set the sample label
#sample_label <- "25164FL-06-01-01_S347"
sample_label <- args[2]

# estimate the library complexity
#estimateLibComplexity(readsDupFreq(bam_file))

# generate and plot the fragment size distribution
jpeg(paste("fragSizeDist", paste(sample_label, "jpg", sep="."), sep="_"))
fragSize <- fragSizeDist(bam_file, sample_label)
dev.off()
