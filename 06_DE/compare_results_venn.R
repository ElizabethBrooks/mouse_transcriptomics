#!/usr/bin/env Rscript

# R script to create venn diagrams comparing results

# import libraries
library(ggVennDiagram)

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/compare_results_venn"
dir.create(workingDir)
setwd(workingDir)

# retrieve previous results
SARA_MANF_BAT_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_SARA_MANF_BAT_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
SARA_MANF_WAT_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_SARA_MANF_WAT_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
SARA_MANF_HYP_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_SARA_MANF_HYP_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
SARA_MANF_LIV_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_SARA_MANF_LIV_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
SARA_BAT_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_SARA_BAT_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
SARA_WAT_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_SARA_WAT_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
SARA_HYP_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_SARA_HYP_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
SARA_LIV_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_SARA_LIV_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
MANF_BAT_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_MANF_BAT_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
MANF_WAT_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_MANF_WAT_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
MANF_HYP_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_MANF_HYP_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")
MANF_LIV_DEGs <- read.csv("/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/treatment_MANF_LIV_M_dge_sig_FDR_0.05_LFC_0.263034405833794.csv")

# create combined lists of DE gene names
SARA_MANF_venn <- list(BAT = SARA_MANF_BAT_DEGs$gene, WAT = SARA_MANF_WAT_DEGs$gene, HYP = SARA_MANF_HYP_DEGs$gene, LIV = SARA_MANF_LIV_DEGs$gene)
SARA_venn <- list(BAT = SARA_BAT_DEGs$gene, WAT = SARA_WAT_DEGs$gene, HYP = SARA_HYP_DEGs$gene, LIV = SARA_LIV_DEGs$gene)
MANF_venn <- list(BAT = MANF_BAT_DEGs$gene, WAT = MANF_WAT_DEGs$gene, HYP = MANF_HYP_DEGs$gene, LIV = MANF_LIV_DEGs$gene)
SARA_MANF_BAT_venn <- list(SARA_MANF = SARA_MANF_BAT_DEGs$gene, SARA = SARA_BAT_DEGs$gene, MANF = MANF_BAT_DEGs$gene)
SARA_MANF_WAT_venn <- list(SARA_MANF = SARA_MANF_WAT_DEGs$gene, SARA = SARA_WAT_DEGs$gene, MANF = MANF_WAT_DEGs$gene)
SARA_MANF_HYP_venn <- list(SARA_MANF = SARA_MANF_HYP_DEGs$gene, SARA = SARA_HYP_DEGs$gene, MANF = MANF_HYP_DEGs$gene)
SARA_MANF_LIV_venn <- list(SARA_MANF = SARA_MANF_LIV_DEGs$gene, SARA = SARA_LIV_DEGs$gene, MANF = MANF_LIV_DEGs$gene)
BAT_venn <- list(SARA = SARA_BAT_DEGs$gene, MANF = MANF_BAT_DEGs$gene)
WAT_venn <- list(SARA = SARA_WAT_DEGs$gene, MANF = MANF_WAT_DEGs$gene)
HYP_venn <- list(SARA = SARA_HYP_DEGs$gene, MANF = MANF_HYP_DEGs$gene)
LIV_venn <- list(SARA = SARA_LIV_DEGs$gene, MANF = MANF_LIV_DEGs$gene)

# create tissue treatment venn diagrams
jpeg("SARA_MANF_results_venn.jpg")
ggVennDiagram(SARA_MANF_venn, label_alpha=0.25, category.names = c("BAT", "WAT", "HYP", "LIV"))
dev.off()
jpeg("SARA_results_venn.jpg")
ggVennDiagram(SARA_venn, label_alpha=0.25, category.names = c("BAT", "WAT", "HYP", "LIV"))
dev.off()
jpeg("MANF_results_venn.jpg")
ggVennDiagram(MANF_venn, label_alpha=0.25, category.names = c("BAT", "WAT", "HYP", "LIV"))
dev.off()
jpeg("SARA_MANF_BAT_results_venn.jpg")
ggVennDiagram(SARA_MANF_BAT_venn, label_alpha=0.25, category.names = c("SARA_MANF", "SARA", "MANF"))
dev.off()
jpeg("SARA_MANF_WAT_results_venn.jpg")
ggVennDiagram(SARA_MANF_WAT_venn, label_alpha=0.25, category.names = c("SARA_MANF", "SARA", "MANF"))
dev.off()
jpeg("SARA_MANF_HYP_results_venn.jpg")
ggVennDiagram(SARA_MANF_HYP_venn, label_alpha=0.25, category.names = c("SARA_MANF", "SARA", "MANF"))
dev.off()
jpeg("SARA_MANF_LIV_results_venn.jpg")
ggVennDiagram(SARA_MANF_LIV_venn, label_alpha=0.25, category.names = c("SARA_MANF", "SARA", "MANF"))
dev.off()
jpeg("BAT_results_venn.jpg")
ggVennDiagram(BAT_venn, label_alpha=0.25, category.names = c("SARA", "MANF"))
dev.off()
jpeg("WAT_results_venn.jpg")
ggVennDiagram(WAT_venn, label_alpha=0.25, category.names = c("SARA", "MANF"))
dev.off()
jpeg("HYP_results_venn.jpg")
ggVennDiagram(HYP_venn, label_alpha=0.25, category.names = c("SARA", "MANF"))
dev.off()
jpeg("LIV_results_venn.jpg")
ggVennDiagram(LIV_venn, label_alpha=0.25, category.names = c("SARA", "MANF"))
dev.off()
