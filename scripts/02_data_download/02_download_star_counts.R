# TCGA Pan-Cancer RNA-seq Download (STAR Counts)
# Project: miR21-SH3GLB1 Prognostic Analysis
# Author: Bhawna Saini
# Date: 2026-02-26


rm(list = ls())
cat("Starting STAR RNA-seq download pipeline...\n")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

cohort_file <- "results/tables/cohort_summary_STAR.csv"
if (!file.exists(cohort_file)) {
  stop("Cohort summary file not found. Run cohort definition first.")
}

cohorts <- read.csv(cohort_file, stringsAsFactors = FALSE)
projects <- cohorts$project

cat("Projects to download:", length(projects), "\n")

for (proj_id in projects) {

   cat("====================================\n")
  cat("Project:", proj_id, "\n")

  tryCatch({

    query <- GDCquery(
      project       = proj_id,
      data.category = "Transcriptome Profiling",
      data.type     = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      access        = "open"
    )

    results <- getResults(query)

    if (nrow(results) == 0) {
      cat("No data found:", proj_id, "\n")
      next
    }

    cat("Files:", nrow(results), "\n")

    GDCdownload(
      query           = query,
      method          = "api",
      directory       = "data/raw",
      files.per.chunk = 6
    )

    cat("Done:", proj_id, "\n")

  }, error = function(e) {
    cat("Failed:", proj_id, "-", conditionMessage(e), "\n")
    write(paste0(Sys.time(), " | ", proj_id, " | ", conditionMessage(e)),
          file = "logs/download_errors.txt", append = TRUE)
  })
}

sink("logs/sessionInfo_download_star.txt")
sessionInfo()
sink()

cat("Finished.\n")
