
###########################################################
# 4) Pathway Enrichment (TCGA)
###########################################################
message("\n=== Running pathway enrichment analysis ===")

if (exists("sig") && nrow(sig) > 0) {
  tryCatch({
    # Prepare gene list
    gene_symbols <- na.omit(sig$symbol)
    
    message(sprintf("Running GO enrichment for %d genes...", length(gene_symbols)))
    
    # GO enrichment
    ego <- enrichGO(
      gene = gene_symbols,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "fdr",
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      write.csv(as.data.frame(ego), "TCGA_GO_enrichment.csv", row.names = FALSE)
      
      pdf("TCGA_GO_dotplot.pdf", width = 12, height = 8)
      print(dotplot(ego, showCategory = 20, title = "GO Enrichment - TCGA-LUAD"))
      dev.off()
      
      pdf("TCGA_GO_barplot.pdf", width = 10, height = 7)
      print(barplot(ego, showCategory = 15, title = "GO Enrichment - TCGA-LUAD"))
      dev.off()
      
      message("Pathway enrichment completed successfully!")
    } else {
      message("No significant GO terms found")
    }
    
  }, error = function(e) {
    message("Error in pathway enrichment: ", e$message)
  })
} else {
  message("No significant genes found for pathway enrichment")
}