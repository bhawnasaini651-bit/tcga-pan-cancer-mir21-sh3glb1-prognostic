sink("logs/01_define_cohorts_STAR.log", split = TRUE)# 01_define_cohorts_STAR.R
# Purpose:
# Identify TCGA cohorts eligible for pan-cancer analysis
# Criteria:
#   - Tumor samples >= 50
#   - Normal samples >= 10
#   - Survival data >= 50 patients
#   - STAR - Counts available

library(TCGAbiolinks)
library(dplyr)

# Get all TCGA projects
projects <- TCGAbiolinks:::getGDCprojects()$project_id
tcga_projects <- projects[grep("TCGA", projects)]

cohort_summary <- data.frame()

for (proj in tcga_projects) {
  
  message("Checking ", proj)
  
  # Query STAR counts
  query <- tryCatch({
    GDCquery(
      project = proj,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
  }, error = function(e) return(NULL))
  
  if (is.null(query)) next
  
  metadata <- getResults(query)
  if (nrow(metadata) == 0) next
  
  sample_type <- metadata$shortLetterCode
  
  tumor_n  <- sum(sample_type == "TP")
  normal_n <- sum(sample_type == "NT")
  
  # Clinical survival data
  clinical <- tryCatch({
    GDCquery_clinic(project = proj, type = "clinical")
  }, error = function(e) return(NULL))
  
  if (is.null(clinical)) next
  
  survival_n <- sum(!is.na(clinical$days_to_death) |
                    !is.na(clinical$days_to_last_follow_up))
  
  cohort_summary <- rbind(
    cohort_summary,
    data.frame(
      project = proj,
      tumor_samples = tumor_n,
      normal_samples = normal_n,
      survival_patients = survival_n
    )
  )
}

# Save full summary
write.csv(cohort_summary,
          "results/tables/cohort_sample_summary.csv",
          row.names = FALSE)

# Apply eligibility filter
eligible_cohorts <- cohort_summary %>%
  filter(
    tumor_samples >= 50,
    normal_samples >= 10,
    survival_patients >= 50
  )

write.csv(eligible_cohorts,
          "results/tables/eligible_cohorts.csv",
          row.names = FALSE)

message("Cohort screening complete.")

sink()
