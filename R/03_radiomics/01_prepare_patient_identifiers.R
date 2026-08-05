# 01_prepare_patient_identifiers.R
#
# Prepare patient identifiers for CT, segmentation, and transcriptomic matching.
#
# Run from an analysis workspace configured through environment variables.

options(stringsAsFactors = FALSE)


# Analysis paths are supplied through environment variables.
required_env_path <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    stop("Environment variable ", name, " is not set.")
  }
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

project_dir <- required_env_path("RADIOGENOMICS_PROJECT_DIR")
setwd(project_dir)

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Missing R package 'openxlsx'. Install dependencies with environment/install_r_dependencies.R.")
}

library(openxlsx)

dirs <- c(
  "02_metadata",
  "05_molecular_scores",
  "07_results"
)

for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

############################################################
# 1. Load molecular score table
############################################################

score_path <- "05_molecular_scores/NSCLC_Radiogenomics_patient_molecular_scores.rds"

if (!file.exists(score_path)) {
  stop(
    "Cannot find molecular score table:\n",
    score_path,
    "\nRun the patient molecular-score step first."
  )
}

score_df <- readRDS(score_path)

cat("\n===== Molecular score table loaded =====\n")
print(dim(score_df))

cat("\n===== Molecular score table columns =====\n")
print(colnames(score_df))

cat("\n===== Sample ID preview =====\n")
print(head(score_df$sampleid, 30))

############################################################
# 2. Standardize patient IDs
############################################################

standardize_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  return(x)
}

score_df$patient_id_molecular <- standardize_id(score_df$sampleid)

# Extract numeric part, useful for checking possible naming variants
score_df$patient_number <- sub("^R01-", "", score_df$patient_id_molecular)

id_table <- data.frame(
  sampleid = score_df$sampleid,
  patient_id_molecular = score_df$patient_id_molecular,
  patient_number = score_df$patient_number,
  CRS = score_df$CRS,
  CRS_z = score_df$CRS_z,
  RSI = score_df$RSI,
  RSI_z = score_df$RSI_z,
  Hypoxia_core10_score = score_df$Hypoxia_core10_score,
  Hypoxia_core10_z = score_df$Hypoxia_core10_z,
  CRS_group_median = score_df$CRS_group_median,
  RSI_group_median = score_df$RSI_group_median,
  Hypoxia_group_median = score_df$Hypoxia_group_median,
  stringsAsFactors = FALSE
)

############################################################
# 3. ID audit
############################################################

id_audit <- data.frame(
  item = c(
    "molecular_patients_total",
    "unique_molecular_patient_ids",
    "duplicated_molecular_patient_ids",
    "example_first_id",
    "example_last_id"
  ),
  value = c(
    nrow(id_table),
    length(unique(id_table$patient_id_molecular)),
    sum(duplicated(id_table$patient_id_molecular)),
    id_table$patient_id_molecular[1],
    id_table$patient_id_molecular[nrow(id_table)]
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Molecular patient ID audit =====\n")
print(id_audit)

cat("\n===== Prepared ID table preview =====\n")
print(head(id_table, 30))

############################################################
# 4. Save outputs
############################################################

write.csv(
  id_table,
  file = "02_metadata/nsclc_radiogenomics_molecular_patient_ids.csv",
  row.names = FALSE
)

saveRDS(
  id_table,
  file = "02_metadata/nsclc_radiogenomics_molecular_patient_ids.rds"
)

wb <- createWorkbook()

addWorksheet(wb, "molecular_patient_ID_table")
writeData(wb, "molecular_patient_ID_table", id_table)

addWorksheet(wb, "ID_audit")
writeData(wb, "ID_audit", id_audit)

saveWorkbook(
  wb,
  file = "02_metadata/nsclc_radiogenomics_molecular_patient_ids.xlsx",
  overwrite = TRUE
)

cat("\n===== Molecular patient-ID table prepared =====\n")
cat("Main output files:\n")
cat("1. 02_metadata/nsclc_radiogenomics_molecular_patient_ids.xlsx\n")
cat("2. 02_metadata/nsclc_radiogenomics_molecular_patient_ids.csv\n")
cat("3. 02_metadata/nsclc_radiogenomics_molecular_patient_ids.rds\n")

cat("1. Molecular patient ID audit\n")
cat("2. Prepared ID table preview\n")
