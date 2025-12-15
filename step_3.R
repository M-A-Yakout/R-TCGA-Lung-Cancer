
###########################################################
# 3) GEO GSE197364 Analysis - FIXED
###########################################################
message("\n=== Starting GEO GSE197364 analysis ===")

tryCatch({
  # Get GEO data
  gse2 <- getGEO("GSE197364", GSEMatrix = TRUE)[[1]]
  expr2 <- exprs(gse2)
  meta2 <- pData(gse2)
  
  # Get group information and clean names
  group_col <- "title"  # Using title column
  
  grp2_raw <- as.character(meta2[[group_col]])
  
  # Simplify group names
  grp2_clean <- gsub("HRSV infected Hep2 20 hpi rep [0-9]+", "HRSV_infected", grp2_raw)
  grp2_clean <- gsub("HRSV infected Hep2 20 hpi rep [0-9]+ \\(UV\\)", "HRSV_UV", grp2_clean)
  grp2_clean <- gsub("Mock infected Hep2 rep [0-9]+", "Mock_control", grp2_clean)
  grp2_clean <- make.names(grp2_clean)
  
  grp2 <- factor(grp2_clean)
  message(sprintf("Groups found: %s", paste(unique(grp2), collapse = ", ")))
  
  # Design matrix
  design2 <- model.matrix(~ 0 + grp2)
  colnames(design2) <- levels(grp2)
  
  # Fit linear model
  fit2 <- lmFit(expr2, design2)
  
  # Make contrasts - HRSV vs Mock
  unique_groups <- levels(grp2)
  if (length(unique_groups) >= 2) {
    # Find HRSV and Mock groups
    hrsv_group <- grep("HRSV_infected", unique_groups, value = TRUE)[1]
    mock_group <- grep("Mock", unique_groups, value = TRUE)[1]
    
    if (!is.na(hrsv_group) && !is.na(mock_group)) {
      cont_formula <- paste(hrsv_group, "-", mock_group)
      message(sprintf("Comparing: %s", cont_formula))
      
      cont_matrix2 <- makeContrasts(
        contrasts = cont_formula,
        levels = design2
      )
      
      fit2 <- contrasts.fit(fit2, cont_matrix2)
      fit2 <- eBayes(fit2)
      
      # Get results
      res2 <- topTable(fit2, number = Inf, adjust.method = "fdr")
      sig2 <- res2 %>% 
        filter(adj.P.Val < 0.05, abs(logFC) >= 1) %>%
        arrange(adj.P.Val)
      
      message(sprintf("Significant DEGs found: %d", nrow(sig2)))
      write.csv(sig2, "GSE197364_sigDEGs.csv", row.names = TRUE)
      
      # Volcano plot
      pdf("GSE197364_volcano.pdf", width = 10, height = 7)
      res2_plot <- res2 %>%
        mutate(
          significant = ifelse(adj.P.Val < 0.05 & abs(logFC) >= 1, "Significant", "Not Significant"),
          neglog10p = -log10(adj.P.Val)
        )
      
      p <- ggplot(res2_plot, aes(x = logFC, y = neglog10p, color = significant)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
        labs(title = sprintf("GSE197364: %s", cont_formula),
             x = "Log Fold Change",
             y = "-Log10 Adjusted P-value") +
        theme_bw() +
        theme(legend.position = "top")
      print(p)
      dev.off()
      
      message("GSE197364 analysis completed successfully!")
    }
  }
  
}, error = function(e) {
  message("Error in GSE197364 analysis: ", e$message)
})