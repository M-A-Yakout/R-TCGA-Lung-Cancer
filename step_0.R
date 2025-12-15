
# 0) Setup and Install Packages
############################################################
message("Installing required packages...")

# Install base package managers
pkgs <- c("BiocManager", "remotes")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

# Bioconductor packages - including apeglm
bioPkgs <- c("TCGAbiolinks", "GEOquery", "DESeq2", "apeglm", "limma", 
             "SummarizedExperiment", "clusterProfiler", "org.Hs.eg.db",
             "pheatmap", "ComplexHeatmap")

need <- bioPkgs[!sapply(bioPkgs, requireNamespace, quietly = TRUE)]
if (length(need) > 0) {
  BiocManager::install(need, ask = FALSE, update = FALSE)
}

# Load libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(GEOquery)
  library(DESeq2)
  library(apeglm)
  library(limma)
  library(SummarizedExperiment)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(dplyr)
  library(ggplot2)
})

# Create results directory
if (!dir.exists("Results")) dir.create("Results")
setwd("Results")