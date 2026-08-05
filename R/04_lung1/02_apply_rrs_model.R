# 02_apply_rrs_model.R
#
# Apply the frozen primary RRS model to Lung1.
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

lung1_root <- required_env_path("LUNG1_DATA_DIR")

############################################################
# 1. Packages
############################################################

required_packages <- c("glmnet", "openxlsx")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install dependencies with environment/install_r_dependencies.R."
  )
}
invisible(lapply(required_packages, library, character.only = TRUE))

library(glmnet)
library(openxlsx)

############################################################
# 2. Paths
############################################################


lung1_feature_csv <- file.path(
  lung1_root,
  "radiomics_features",
  "lung1_radiomics_features.csv"
)

lung1_feature_log_csv <- file.path(
  lung1_root,
  "radiomics_features",
  "lung1_radiomics_features_log.csv"
)

clinical_rds <- file.path(
  lung1_root,
  "clinical",
  "Lung1_clinical_clean_initial.rds"
)

frozen_model_rds <- file.path(
  project_dir,
  "07_results",
  "rrs_model",
  "rrs_elastic_net_model.rds"
)

out_dir <- file.path(
  lung1_root,
  "external_RRS"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(lung1_feature_csv)) {
  stop("Cannot find Lung1 feature CSV: ", lung1_feature_csv)
}

if (!file.exists(lung1_feature_log_csv)) {
  stop("Cannot find Lung1 feature log CSV: ", lung1_feature_log_csv)
}

if (!file.exists(clinical_rds)) {
  stop("Cannot find Lung1 clinical RDS: ", clinical_rds)
}

if (!file.exists(frozen_model_rds)) {
  stop("Cannot find the fitted RRS model object: ", frozen_model_rds)
}

############################################################
# 3. Load data
############################################################

cat("\n===== Loading Lung1 PyRadiomics features =====\n")

feat <- read.csv(
  lung1_feature_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

feat_log <- read.csv(
  lung1_feature_log_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Feature table dim:\n")
print(dim(feat))

cat("\nFeature extraction status:\n")
print(table(feat_log$status, useNA = "ifany"))

cat("\nFailed feature extraction patients:\n")
print(feat_log[feat_log$status == "FAILED", c("patient_id", "message")])

cat("\n===== Loading clinical table =====\n")

clinical <- readRDS(clinical_rds)

cat("Clinical table dim:\n")
print(dim(clinical))

cat("\nClinical columns:\n")
print(colnames(clinical))

cat("\n===== Loading the fitted RRS model =====\n")

frozen_model <- readRDS(frozen_model_rds)

cat("Frozen model names:\n")
print(names(frozen_model))

cat("\nFrozen model final alpha:\n")
print(frozen_model$final_alpha)

cat("\nFrozen model final lambda:\n")
print(frozen_model$final_lambda)

cat("\nFrozen model feature number:\n")
print(length(frozen_model$final_features))

############################################################
# 4. Standardize IDs
############################################################

feat$patient_id <- toupper(trimws(as.character(feat$patient_id)))
feat_log$patient_id <- toupper(trimws(as.character(feat_log$patient_id)))

if (!"patient_id" %in% colnames(clinical)) {
  if ("PatientID" %in% colnames(clinical)) {
    clinical$patient_id <- toupper(trimws(as.character(clinical$PatientID)))
  } else {
    stop("Cannot find patient_id or PatientID in clinical table.")
  }
}

clinical$patient_id <- toupper(trimws(as.character(clinical$patient_id)))

############################################################
# 5. Feature audit
############################################################

radiomics_cols <- grep("^original_", colnames(feat), value = TRUE)

cat("\n===== Lung1 radiomics feature audit =====\n")
cat("Original radiomics features in Lung1:\n")
print(length(radiomics_cols))

cat("\nFrozen model required features:\n")
print(length(frozen_model$final_features))

missing_model_features <- setdiff(frozen_model$final_features, colnames(feat))
extra_lung1_features <- setdiff(radiomics_cols, frozen_model$final_features)

cat("\nMissing frozen model features in Lung1:\n")
print(missing_model_features)

cat("\nExtra Lung1 original features not used by frozen model:\n")
print(length(extra_lung1_features))

if (length(missing_model_features) > 0) {
  stop(
    "Some frozen model features are missing in Lung1. Cannot apply model. Missing: ",
    paste(missing_model_features, collapse = "; ")
  )
}

############################################################
# 6. Prepare model matrix using frozen features only
############################################################

model_features <- frozen_model$final_features

x_lung1 <- feat[, model_features, drop = FALSE]

for (cc in model_features) {
  x_lung1[[cc]] <- suppressWarnings(as.numeric(x_lung1[[cc]]))
}

feature_missing_summary <- data.frame(
  feature = model_features,
  missing_n = sapply(model_features, function(cc) sum(is.na(x_lung1[[cc]]))),
  infinite_n = sapply(model_features, function(cc) sum(is.infinite(x_lung1[[cc]]))),
  stringsAsFactors = FALSE
)

cat("\n===== External feature missingness for frozen model features =====\n")
print(feature_missing_summary)

if (any(feature_missing_summary$missing_n > 0) || any(feature_missing_summary$infinite_n > 0)) {
  stop("Missing or infinite values detected in Lung1 frozen model features.")
}

x_lung1_mat <- as.matrix(x_lung1)

############################################################
# 7. Apply frozen center and scale from training cohort
############################################################

center <- frozen_model$center
scale <- frozen_model$scale

if (!all(model_features %in% names(center))) {
  stop("Frozen center vector does not contain all model features.")
}

if (!all(model_features %in% names(scale))) {
  stop("Frozen scale vector does not contain all model features.")
}

center <- center[model_features]
scale <- scale[model_features]

scale[is.na(scale) | scale == 0] <- 1

x_lung1_scaled <- sweep(x_lung1_mat, 2, center, "-")
x_lung1_scaled <- sweep(x_lung1_scaled, 2, scale, "/")

############################################################
# 8. Predict external RRS
############################################################

rrs_pred <- as.numeric(
  predict(
    frozen_model$final_cvfit,
    newx = x_lung1_scaled,
    s = frozen_model$final_lambda
  )
)

rrs_df <- data.frame(
  patient_id = feat$patient_id,
  Lung1_RRS = rrs_pred,
  stringsAsFactors = FALSE
)

rrs_df$Lung1_RRS_z <- as.numeric(scale(rrs_df$Lung1_RRS))
rrs_df$Lung1_RRS_group_median <- ifelse(
  rrs_df$Lung1_RRS >= median(rrs_df$Lung1_RRS, na.rm = TRUE),
  "RRS_high",
  "RRS_low"
)

cat("\n===== Lung1 external RRS summary =====\n")
print(summary(rrs_df$Lung1_RRS))

cat("\nRRS group table:\n")
print(table(rrs_df$Lung1_RRS_group_median, useNA = "ifany"))

############################################################
# 9. Merge RRS with clinical OS table
############################################################

analysis_df <- merge(
  clinical,
  rrs_df,
  by = "patient_id",
  all = FALSE
)

cat("\n===== Lung1 RRS + clinical merged table =====\n")
print(dim(analysis_df))

cat("\nMerged columns:\n")
print(colnames(analysis_df))

############################################################
# 10. Standardize survival columns
############################################################

if (!"Survival.time" %in% colnames(analysis_df)) {
  stop("Survival.time column not found.")
}

if (!"deadstatus.event" %in% colnames(analysis_df)) {
  stop("deadstatus.event column not found.")
}

analysis_df$OS_time_days <- suppressWarnings(as.numeric(analysis_df$Survival.time))
analysis_df$OS_event <- suppressWarnings(as.numeric(analysis_df$deadstatus.event))

analysis_df$OS_time_months <- analysis_df$OS_time_days / 30.4375

analysis_df <- analysis_df[
  is.finite(analysis_df$OS_time_days) &
    is.finite(analysis_df$OS_event) &
    analysis_df$OS_time_days > 0 &
    analysis_df$OS_event %in% c(0, 1),
]

cat("\n===== Survival variable audit =====\n")
cat("Analysis patients after OS filtering:\n")
print(nrow(analysis_df))

cat("\nOS event table:\n")
print(table(analysis_df$OS_event, useNA = "ifany"))

cat("\nOS time days summary:\n")
print(summary(analysis_df$OS_time_days))

cat("\nOS time months summary:\n")
print(summary(analysis_df$OS_time_months))

############################################################
# 11. Sample flow summary
############################################################

sample_flow <- data.frame(
  step = c(
    "Lung1 clinical patients",
    "Downloaded CT+SEG patients",
    "Converted final tumor-only masks PASS",
    "PyRadiomics extraction SUCCESS",
    "External RRS calculated",
    "RRS + OS analysis patients",
    "Excluded during PyRadiomics"
  ),
  n = c(
    length(unique(clinical$patient_id)),
    421,
    421,
    nrow(feat),
    nrow(rrs_df),
    nrow(analysis_df),
    sum(feat_log$status == "FAILED", na.rm = TRUE)
  ),
  note = c(
    "Original Lung1 clinical table",
    "LUNG1-128 missing CT+SEG download",
    "Final tumor-only mask QC PASS",
    "3 cases failed due to single-voxel tumor mask",
    "Fitted RRS model applied",
    "Complete OS and RRS",
    paste(feat_log$patient_id[feat_log$status == "FAILED"], collapse = "; ")
  ),
  stringsAsFactors = FALSE
)

cat("\n===== Sample flow =====\n")
print(sample_flow)

############################################################
# 12. Save outputs
############################################################

write.csv(
  rrs_df,
  file = file.path(out_dir, "lung1_rrs_scores.csv"),
  row.names = FALSE
)

write.csv(
  analysis_df,
  file = file.path(out_dir, "lung1_rrs_survival_dataset.csv"),
  row.names = FALSE
)

write.csv(
  sample_flow,
  file = file.path(out_dir, "lung1_rrs_sample_flow.csv"),
  row.names = FALSE
)

write.csv(
  feature_missing_summary,
  file = file.path(out_dir, "lung1_rrs_feature_missingness.csv"),
  row.names = FALSE
)

saveRDS(
  rrs_df,
  file = file.path(out_dir, "lung1_rrs_scores.rds")
)

saveRDS(
  analysis_df,
  file = file.path(out_dir, "lung1_rrs_survival_dataset.rds")
)

wb <- createWorkbook()

addWorksheet(wb, "sample_flow")
writeData(wb, "sample_flow", sample_flow)

addWorksheet(wb, "RRS_summary")
writeData(wb, "RRS_summary", rrs_df)

addWorksheet(wb, "analysis_table")
writeData(wb, "analysis_table", analysis_df)

addWorksheet(wb, "feature_missingness")
writeData(wb, "feature_missingness", feature_missing_summary)

addWorksheet(wb, "feature_failed_cases")
writeData(wb, "feature_failed_cases", feat_log[feat_log$status == "FAILED", ])

saveWorkbook(
  wb,
  file = file.path(out_dir, "lung1_rrs_application_summary.xlsx"),
  overwrite = TRUE
)

cat("\n===== Lung1 RRS application completed =====\n")
cat("Main output files:\n")

cat("1. Lung1 external RRS summary\n")
cat("2. RRS group table\n")
cat("3. Lung1 RRS + clinical merged table\n")
cat("4. Survival variable audit\n")
cat("5. Sample flow\n")
