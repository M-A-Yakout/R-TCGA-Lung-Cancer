
###########################################################
# Summary
###########################################################
message("\n", paste(rep("=", 60), collapse = ""))
message("Analysis Complete!")
message("Results saved in: ", getwd())
message(paste(rep("=", 60), collapse = ""), "\n")

# List all generated files
files <- list.files(pattern = "\\.(csv|pdf|rds)$")
if (length(files) > 0) {
  message("Generated files:")
  for (f in files) message("  - ", f)
} else {
  message("No output files were generated. Check errors above.")
}