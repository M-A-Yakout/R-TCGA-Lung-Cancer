
###########################################################
# 1) TCGA-LUAD Complete Analysis (STAR counts)
###########################################################
message("\n=== Starting TCGA-LUAD analysis ===")

tryCatch({
  # Query TCGA data
  query <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
  )
  
  # Download data
  message("Downloading TCGA data...")
  GDCdownload(query, method = "api", files.per.chunk = 10)
  
  # Prepare data
  message("Preparing TCGA data...")
  se <- GDCprepare(query, summarizedExperiment = TRUE)
  cnt <- assay(se, "unstranded")
  
  # Define groups based on sample type
  sample_types <- colData(se)$sample_type
  grp <- factor(
    ifelse(sample_types == "Primary Tumor", "Tumor", "Normal"),
    levels = c("Normal", "Tumor")
  )
  
  # Filter samples
  keep <- !is.na(grp)
  cnt <- cnt[, keep]
  grp <- grp[keep]
  
  message(sprintf("Samples: %d Normal, %d Tumor", 
                  sum(grp == "Normal"), sum(grp == "Tumor")))
  
  # Create DESeq2 dataset
  coldata <- DataFrame(condition = grp, row.names = colnames(cnt))
  
  # Filter low count genes
  cnt <- cnt[rowSums(cnt >= 10) >= 3, ]
  message(sprintf("Genes after filtering: %d", nrow(cnt)))
  
  # DESeq2 analysis
  message("Running DESeq2 analysis...")
  dds <- DESeqDataSetFromMatrix(cnt, coldata, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
  
  # Shrink log fold changes with apeglm
  message("Shrinking log fold changes...")
  res <- lfcShrink(dds, coef = 2, res = res, type = "apeglm")
  
  # Extract significant DEGs
  sig <- as.data.frame(res) %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) %>%
    mutate(
      ensembl = sub("\\..*", "", rownames(.)),
      symbol = mapIds(org.Hs.eg.db, ensembl, "SYMBOL", "ENSEMBL", multiVals = "first")
    ) %>%
    arrange(padj)
  
  message(sprintf("Significant DEGs found: %d", nrow(sig)))
  write.csv(sig, "TCGA_LUAD_sigDEGs.csv", row.names = FALSE)
  
  # Custom volcano plot
  message("Creating volcano plot...")
  res_df <- as.data.frame(res) %>%
    mutate(
      significant = ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1, 
                           "Significant", "Not Significant"),
      neglog10padj = -log10(padj)
    )
  
  pdf("TCGA_LUAD_volcano.pdf", width = 10, height = 7)
  p <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10padj, color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    labs(title = "TCGA-LUAD: Tumor vs Normal",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_bw() +
    theme(legend.position = "top")
  print(p)
  dev.off()
  
  message("TCGA-LUAD volcano plot saved!")
  
  # Save results object for later use
  saveRDS(res, "TCGA_LUAD_results.rds")
  saveRDS(dds, "TCGA_LUAD_dds.rds")
  
  message("TCGA-LUAD analysis completed successfully!")
  
}, error = function(e) {
  message("Error in TCGA analysis: ", e$message)
})