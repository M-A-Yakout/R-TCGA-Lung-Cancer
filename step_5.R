
###########################################################
# 5) Heatmap - Top 50 DEGs (TCGA)
###########################################################
message("\n=== Creating heatmap for top 50 DEGs ===")

if (exists("sig") && exists("dds") && nrow(sig) >= 50) {
  tryCatch({
    # Get top 50 genes
    top50_genes <- head(rownames(sig), 50)
    
    # Get normalized counts
    norm_counts <- counts(dds, normalized = TRUE)
    top50_data <- norm_counts[top50_genes, ]
    
    # Add gene symbols to row names
    gene_names <- sig$symbol[match(rownames(top50_data), rownames(sig))]
    rownames(top50_data) <- make.unique(as.character(gene_names))
    
    # Annotation
    ann_col <- data.frame(
      Group = grp,
      row.names = colnames(top50_data)
    )
    
    # Create heatmap
    pdf("TCGA_LUAD_heatmap_top50.pdf", width = 12, height = 14)
    pheatmap(
      log2(top50_data + 1),
      scale = "row",
      annotation_col = ann_col,
      show_colnames = FALSE,
      cluster_cols = TRUE,
      cluster_rows = TRUE,
      main = "Top 50 DEGs - TCGA-LUAD (Tumor vs Normal)",
      fontsize_row = 7,
      fontsize = 10,
      color = colorRampPalette(c("blue", "white", "red"))(100)
    )
    dev.off()
    
    message("Heatmap created successfully!")
    
  }, error = function(e) {
    message("Error creating heatmap: ", e$message)
  })
} else if (exists("sig") && nrow(sig) > 0 && nrow(sig) < 50) {
  message(sprintf("Only %d significant genes found (need 50 for heatmap)", nrow(sig)))
} else {
  message("No significant genes found for heatmap")
}
