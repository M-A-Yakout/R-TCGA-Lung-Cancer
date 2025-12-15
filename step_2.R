
###########################################################
# 2) GEO GSE188427 Analysis - FIXED
###########################################################
message("\n=== Starting GEO GSE188427 analysis ===")

tryCatch({
  # Get GEO data
  gse1 <- getGEO("GSE188427", GSEMatrix = TRUE)[[1]]
  expr1 <- exprs(gse1)
  meta1 <- pData(gse1)
  
  # Get group information and clean names
  group_col <- "characteristics_ch1.2"  # Adjust if needed
  
  if (!group_col %in% colnames(meta1)) {
    message("Available metadata columns:")
    print(head(colnames(meta1), 20))
    group_col <- grep("source_name|characteristics|group|condition", 
                      colnames(meta1), ignore.case = TRUE, value = TRUE)[1]
  }
  
  # Create clean group names
  grp1_raw <- as.character(meta1[[group_col]])
  
  # Simplify group names (e.g., "whole blood_RSV_Day 1_IN" -> "RSV_Day1_IN")
  grp1_clean <- gsub("whole blood_", "", grp1_raw)
  grp1_clean <- gsub("Healthy person", "Healthy", grp1_clean)
  grp1_clean <- gsub(" ", "_", grp1_clean)
  grp1_clean <- make.names(grp1_clean)  # Make R-valid names
  
  grp1 <- factor(grp1_clean)
  message(sprintf("Groups found: %s", paste(unique(grp1), collapse = ", ")))
  
  # Design matrix
  design1 <- model.matrix(~ 0 + grp1)
  colnames(design1) <- levels(grp1)
  
  # Fit linear model
  fit1 <- lmFit(expr1, design1)
  
  # Make contrasts - compare first two unique groups
  unique_groups <- levels(grp1)
  if (length(unique_groups) >= 2) {
    cont_formula <- paste(unique_groups[2], "-", unique_groups[1])
    message(sprintf("Comparing: %s", cont_formula))
    
    cont_matrix1 <- makeContrasts(
      contrasts = cont_formula,
      levels = design1
    )
    
    fit1 <- contrasts.fit(fit1, cont_matrix1)
    fit1 <- eBayes(fit1)
    
    # Get results
    res1 <- topTable(fit1, number = Inf, adjust.method = "fdr")
    sig1 <- res1 %>% 
      filter(adj.P.Val < 0.05, abs(logFC) >= 1) %>%
      arrange(adj.P.Val)
    
    message(sprintf("Significant DEGs found: %d", nrow(sig1)))
    write.csv(sig1, "GSE188427_sigDEGs.csv", row.names = TRUE)
    
    # Volcano plot
    pdf("GSE188427_volcano.pdf", width = 10, height = 7)
    res1_plot <- res1 %>%
      mutate(
        significant = ifelse(adj.P.Val < 0.05 & abs(logFC) >= 1, "Significant", "Not Significant"),
        neglog10p = -log10(adj.P.Val)
      )
    
    p <- ggplot(res1_plot, aes(x = logFC, y = neglog10p, color = significant)) +
      geom_point(alpha = 0.5, size = 1) +
      scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      labs(title = sprintf("GSE188427: %s", cont_formula),
           x = "Log Fold Change",
           y = "-Log10 Adjusted P-value") +
      theme_bw() +
      theme(legend.position = "top")
    print(p)
    dev.off()
    
    message("GSE188427 analysis completed successfully!")
  }
  
}, error = function(e) {
  message("Error in GSE188427 analysis: ", e$message)
})