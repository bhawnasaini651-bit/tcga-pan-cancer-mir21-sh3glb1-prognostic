############################################################
# TCGA-BRCA Expression + Survival Prototype
# Project: miR21–SH3GLB1 Prognostic Analysis
############################################################

rm(list = ls())
cat("Starting BRCA expression pipeline...\n")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(survival)
library(survminer)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

project_id <- "TCGA-BRCA"

cat("Querying BRCA STAR counts...\n")

query <- GDCquery(
  project = project_id,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

cat("Extracting expression matrix...\n")

expr_matrix <- assay(data)

# Extract genes of interest
genes_of_interest <- c("MIR21", "SH3GLB1")

gene_data <- expr_matrix[rownames(expr_matrix) %in% genes_of_interest, ]

gene_df <- as.data.frame(t(gene_data))
gene_df$barcode <- rownames(gene_df)

# Clinical data
clinical <- colData(data) %>% as.data.frame()
clinical$barcode <- rownames(clinical)

merged <- inner_join(gene_df, clinical, by = "barcode")

cat("Preparing survival object...\n")

merged$OS.time <- as.numeric(merged$days_to_death)
merged$OS.time[is.na(merged$OS.time)] <- as.numeric(merged$days_to_last_follow_up[is.na(merged$OS.time)])

merged$OS <- ifelse(merged$vital_status == "Dead", 1, 0)

# Median split on MIR21
merged$MIR21_group <- ifelse(
  merged$MIR21 > median(merged$MIR21, na.rm = TRUE),
  "High",
  "Low"
)

surv_object <- Surv(merged$OS.time, merged$OS)

fit <- survfit(surv_object ~ MIR21_group, data = merged)

plot <- ggsurvplot(
  fit,
  data = merged,
  pval = TRUE,
  risk.table = TRUE,
  title = "TCGA-BRCA: MIR21 Survival"
)

ggsave("results/plots/BRCA_MIR21_survival.png", plot$plot, width = 6, height = 5)

write.csv(merged, "data/processed/BRCA_expression_clinical.csv", row.names = FALSE)

sink("logs/sessionInfo_brca_expression.txt")
sessionInfo()
sink()

cat("BRCA pipeline complete.\n")
