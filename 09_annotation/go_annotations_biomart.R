#!/usr/bin/env Rscript

# script to retrieve mus musculus annotations

# load libraries
library(biomaRt)
library(tibble)

# set working directory
workingDir="/Users/bamflappy/MackLab/metabolic_adaptation/annotations"
dir.create(workingDir)
setwd(workingDir)

# import gene count data
inputTable <- read.csv(file="/Users/bamflappy/MackLab/metabolic_adaptation/counted/counts_merged_samples.csv", row.names="gene")

# import normalized counts
normTable <- read.csv(file="/Users/bamflappy/MackLab/metabolic_adaptation/Biostatistics/DEAnalysis_14Sep2026/single_means_random/normalizedCounts.csv", row.names="gene")

# counts gene list
counted_genes <- rownames(inputTable)
norm_genes <- rownames(normTable)

# Ensembl BioMart database
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = "useast")

# list attributes
#attributes <- listAttributes(ensembl)

# list filters
#filters <- listFilters(ensembl)

# get GO IDs
go_terms <- getBM(
  attributes = c("ensembl_gene_id", "go_id", "gene_biotype"),
  mart = ensembl
)

# export annotations
go_terms_list <- as_tibble(go_terms)
write.table(go_terms_list, file="go_terms_biomaRt.csv", sep=",", row.names=FALSE, quote=FALSE)

# filter annotations to genes in the counts file
go_terms_counted <- go_terms[go_terms$ensembl_gene_id %in% counted_genes,]
go_terms_norm <- go_terms[go_terms$ensembl_gene_id %in% norm_genes,]

# check how many genes do not have annotations
length(counted_genes) - length(unique(go_terms_counted$ensembl_gene_id))
length(norm_genes) - length(unique(go_terms_norm$ensembl_gene_id))

# get the protein coding genes
go_terms_coding <- go_terms[go_terms$gene_biotype == "protein_coding",]

# export annotations
go_terms_coding$gene_biotype <- NULL
go_terms_coding_list <- as_tibble(go_terms_coding)
write.table(go_terms_coding_list, file="go_terms_coding_biomaRt.csv", sep=",", row.names=FALSE, quote=FALSE)

# filter annotations to genes in the counts file
go_terms_coding_counted <- go_terms_counted[go_terms_counted$gene_biotype == "protein_coding",]
go_terms_coding_norm <- go_terms_norm[go_terms_norm$gene_biotype == "protein_coding",]

# check how many coding genes do not have annotations
length(counted_genes) - length(unique(go_terms_coding_counted$ensembl_gene_id))
length(norm_genes) - length(unique(go_terms_coding_norm$ensembl_gene_id))
