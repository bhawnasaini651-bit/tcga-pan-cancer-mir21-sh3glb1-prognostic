#!/usr/bin/env Rscript

# ==========================================================
# TCGA Pan-Cancer Cohort Definition Script (STAR Counts)
# Reproducible + Git-Ready Version
# ==========================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(dplyr)
})

# -----------------------------
# Setup
# -----------------------------
start_time <- Sys.time()

dir.create("logs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

log_file <- "logs/01_define_cohorts_STAR.log"
sink(log_file, split = TRUE)

cat("====================================================\n")
cat("TCGA Cohort Definition (STAR Counts)\n")
cat("Execution time:", format(start_time), "\n")
cat("====================================================\n\n")

# -----------------------------
# Retrieve TCGA Projects
# -----------------------------
cat("Retrieving TCGA project list...\n")

projects_df <- getGDCprojects()

tcga_projects <- projects_df %>%
  filter(grepl("^TCGA-", project_id)) %>%
  pull(project_id) %>%
  sort()

cat("Total TCGA projects detected:", length(tcga_projects), "\n\n")

# -----------------------------
# Initialize container
# -----------------------------
cohort_summary <- data.frame()

# -----------------------------
# Loop Through Projects
# -----------------------------
for (tcga_project in tcga_projects) {

  cat("----------------------------------------------------\n")
  cat("Processing:", tcga_project, "\n")

  query <- tryCatch({
    GDCquery(
      project = tcga_project,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
  }, error = function(e) {
    cat("Query failed:", e$message, "\n")
    return(NULL)
  })

  if (is.null(query)) next

  results <- tryCatch({
    getResults(query)
  }, error = function(e) {
    cat("Failed retrieving results:", e$message, "\n")
    return(NULL)
  })

  if (is.null(results) || nrow(results) == 0) {
    cat("No STAR files found.\n")
    next
  }

  # -----------------------------
  # Sample Counting
  # -----------------------------
  sample_type <- results$sample_type

  tumor_n  <- sum(sample_type == "Primary Tumor", na.rm = TRUE)
  normal_n <- sum(sample_type == "Solid Tissue Normal", na.rm = TRUE)

  cat("Tumor samples:", tumor_n, "\n")
  cat("Normal samples:", normal_n, "\n")
  cat("Total files:", nrow(results), "\n")

  # -----------------------------
  # Filtering Threshold
  # (Adjust based on survival needs)
  # -----------------------------
  if (tumor_n < 50) {
    cat("Skipped (tumor_n < 50)\n")
    next
  }

  # -----------------------------
  # Append to summary
  # -----------------------------
  cohort_summary <- rbind(
    cohort_summary,
    data.frame(
      project = tcga_project,
      tumor_samples = tumor_n,
      normal_samples = normal_n,
      total_files = nrow(results),
      stringsAsFactors = FALSE
    )
  )
}

# -----------------------------
# Final Validation
# -----------------------------
if (nrow(cohort_summary) == 0) {
  cat("\nNo cohorts passed filtering threshold.\n")
  sink()
  stop("No cohorts were successfully processed.")
}

# -----------------------------
# Save Results
# -----------------------------
output_file <- "results/tables/cohort_summary_STAR.csv"

write.csv(
  cohort_summary,
  output_file,
  row.names = FALSE
)

cat("\n====================================================\n")
cat("Cohort definition completed successfully.\n")
cat("Cohorts passing threshold:", nrow(cohort_summary), "\n")
cat("Results saved to:", output_file, "\n")
cat("====================================================\n")

# -----------------------------
# Save Session Info (Reproducibility)
# -----------------------------
session_file <- "logs/sessionInfo_define_cohorts.txt"
writeLines(capture.output(sessionInfo()), session_file)

end_time <- Sys.time()
cat("Execution finished at:", format(end_time), "\n")
cat("Total runtime:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

sink()
