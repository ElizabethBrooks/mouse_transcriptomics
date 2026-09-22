#!/usr/bin/env Rscript

# R script to visualize the counting statistics from HTSeq

# turn off scientific notation
options(scipen = 999)

# import libraries
library(tibble)
library(rcartocolor)
library(dplyr)
library(ggplot2)

# plotting Palettes
# https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
# https://github.com/Nowosad/rcartocolor
plotColors <- carto_pal(12, "Safe")
plotColorSubset <- c(plotColors[4], plotColors[5], plotColors[6])
plotColorSubset <- c(plotColors[3], plotColors[4], plotColors[5], plotColors[6])

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/stats_counted"
dir.create(workingDir)
setwd(workingDir)

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/MackLab/metabolic_adaptation/counted/counts_merged_samples.csv", row.names="gene")

# import grouping factor
targets <- read.csv(file="/Users/bamflappy/MackLab/input_data/experimental_design.csv", row.names="sample")

# retrieve the lines with counting statistics (htseq)
keep <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
countsTable <- inputTable[row.names(inputTable) %in% keep,]

# retrieve read counts for each sample
sample_read_counts <- colSums(inputTable)

# subset the data frame for transformation
countsTable1 <- data.frame(
  sample = colnames(countsTable),
  stat = "no_feature",
  count = countsTable[,1],
  percent = countsTable[,1] / sample_read_counts * 100
)
countsTable2 <- data.frame(
  sample = colnames(countsTable),
  stat = "ambiguous",
  count = countsTable[,2],
  percent = countsTable[,2] / sample_read_counts * 100
)
countsTable3 <- data.frame(
  sample = colnames(countsTable),
  stat = "too_low_aQual",
  count = countsTable[,3],
  percent = countsTable[,3] / sample_read_counts * 100
)
countsTable4 <- data.frame(
  sample = colnames(countsTable),
  stat = "not_aligned",
  count = countsTable[,4],
  percent = countsTable[,4] / sample_read_counts * 100
)
countsTable5 <- data.frame(
  sample = colnames(countsTable),
  stat = "alignment_not_unique",
  count = countsTable[,5],
  percent = countsTable[,5] / sample_read_counts * 100
)

# transform the data for plotting
countsTable_long <- rbind(countsTable1, countsTable2, countsTable3, countsTable4, countsTable5)

# violin plot of stat counts
counts_plot <- ggplot(countsTable_long, aes(x=stat, y=count, fill=stat)) + 
  geom_violin() +
  #xlab("Statistic") +
  ylab("Count") +
  theme(legend.position="none") +
  coord_flip() +
  theme_bw(base_size = 18)
# save the plot to a png file
ggsave("counted_stats.png", plot = counts_plot, bg = "white", device = "png", width = 16, height = 8, units = "in")

# violin plot of stat percents
percent_plot <- ggplot(countsTable_long, aes(x=stat, y=percent, fill=stat)) + 
  geom_violin() +
  xlab("Statistic") +
  ylab("Percent of Quantified Reads Per Sample") +
  theme(legend.position="none") +
  coord_flip() +
  theme_bw(base_size = 18)
# save the plot to a png file
ggsave("counted_stats_percent.png", plot = percent_plot, bg = "white", device = "png", width = 16, height = 8, units = "in")
